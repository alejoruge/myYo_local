function main_process()
% Processes 40x OCT images with optional distortion and coverslip removal.
% USAGE:
%   main_process()
%   Configures input file path and processing options via internal variables.
% INPUTS:
%   tifFilePath - Path to the .tiff file for processing.
%   skipDistortionRemoval - true to skip distortion correction, false to apply it.
%   skipCoverslipRemoval - true to skip coverslip artifact removal, false to apply it.
%   pixel_size_um - Microns per pixel for heatmap scaling.
% OUTPUTS:
%   Outputs are not returned, but files are saved based on processing steps applied.
%   Also, it will always create a heatmap to the surface of the tissue.

% Internally calls:
%   distortion_removal() - Corrects parabolic distortion.
%   process_oct_data() - Detects surface, modifies data, and generates heatmap.
%   Various utility functions handle specific processing steps.

%% User Inputs

    % Set the TIF file path for processing
    tifFilePath = '\';

    % Set to true if you want to skip parabolic distortion removal
    skipDistortionRemoval = false;

    % Set to true if you want to skip coverslip removal
    skipCoverslipRemoval = false;

    % Microns per pixel for scaling the heatmap (only heatmap)
    pixel_size_um = 1;


%% Process the image based on input
    
    if ~skipDistortionRemoval
        [data, metadata, clim] = distortion_removal(tifFilePath);  % Corrected data if not skipped
    else
        [data, metadata, clim] = yOCTFromTif(tifFilePath);  % Load data directly if skipped
    end
    
%% Proceed to modify OCT data

    outputFilePath = process_oct_data(data, metadata, clim, pixel_size_um, tifFilePath, skipCoverslipRemoval, skipDistortionRemoval);

    if ~isempty(outputFilePath)
        fprintf('Final modified file saved at %s\n', outputFilePath);
    else
        fprintf('No modifications made. Only heatmap generated and no new file created.\n');
    end
end


%% Internal functions here:

function [data, metadata, clim] = distortion_removal(tifFilePath)
% Corrects parabolic distortion in OCT images using polynomial coefficients tailored for a 40x lens.

    % Polynomial coefficients for OCT 40x
    pXCenter = [2.9625e-12, -4.1659e-09, 2.0480e-06, -0.0011, 0.3646, 64.6771];
    pYCenter = [2.2165e-14, 3.5006e-11, -4.9389e-08, -4.7830e-04, 0.2484, 76.1154];

    % Load the image data
    [data, metadata, clim] = yOCTFromTif(tifFilePath);
    [zSize, xSize, ySize] = size(data);

    % Apply polynomial correction to each pixel
    [X, Y] = meshgrid(1:xSize, 1:ySize);
    Z_pred = polyval(pXCenter, X) + polyval(pYCenter, Y);
    centerShift = round(Z_pred(floor(ySize/2), floor(xSize/2)));

    correctedData = zeros(size(data), 'like', data);
    for x = 1:xSize
        for y = 1:ySize
            zShift = round(Z_pred(y, x) - centerShift);
            correctedData(:, x, y) = circshift(data(:, x, y), -zShift);
        end
    end
    
    data = correctedData;  % Update data with corrected data
    fprintf('Distortion correction applied.\n');
    return;
end

%% Delete coverslip and detect surface:
function outputFilePath = process_oct_data(data, metadata, clim, pixel_size_um, tifFilePath, skipCoverslipRemoval, skipDistortionRemoval)
    [z_size, x_size, y_size] = size(data);
    disp(['-- Data dimensions: ', num2str(z_size), 'x', num2str(x_size), 'x', num2str(y_size)]);
    
    % Initialize an array to hold surface depth for each (x, y) coordinate
    surface_depth = zeros(x_size, y_size);

    % Process each slice to find the surface depth
    for y = 1:y_size
        current_slice = data(:, :, y);
        surface_depth(:, y) = detect_surface_depth(current_slice, z_size);
    end

    % Verify the surface depth matrix
    disp(['-- Surface depth matrix min: ', num2str(min(surface_depth(:))), ', max: ', num2str(max(surface_depth(:)))]);

    % Apply Gaussian smoothing to the surface depth data
    sigma = 1.5;  % Standard deviation for the Gaussian kernel
    smoothed_surface_depth = round(smooth_surface(surface_depth, sigma));

    % Display the surface depth as a heatmap
    figure;
    imagesc((smoothed_surface_depth * pixel_size_um).');
    colormap(flipud(jet));
    colorbar;
    set(gca, 'YDir', 'normal');
    title('Surface Depth Heatmap');
    xlabel('X-axis');
    ylabel('Y-axis');

    % Modify the data based on the detected surface depth if not skipping
    if ~skipCoverslipRemoval
        modified_data = modify_data_based_on_surface(data, smoothed_surface_depth);
        disp(['-- Modified data min/max: ', num2str(min(modified_data(:))), '/', num2str(max(modified_data(:)))]);
    else
        modified_data = data;  % Skip modification, use original data
        disp('-- Modification step skipped (no coverslip removal).');
    end

    % Verify the modification
    disp(['-- Modified data min: ', num2str(min(modified_data(:))), ', max: ', num2str(max(modified_data(:)))]);

    % Generate file name based on the operations performed
    outputFilePath = '';
    if ~skipDistortionRemoval || ~skipCoverslipRemoval
        [pathstr, name, ~] = fileparts(tifFilePath);
        newFileName = name;
        if ~skipDistortionRemoval && ~skipCoverslipRemoval
            newFileName = [newFileName '_DistCorrected_CoverslipRemoved'];
        elseif ~skipDistortionRemoval
            newFileName = [newFileName '_DistCorrected'];
        elseif ~skipCoverslipRemoval
            newFileName = [newFileName '_CoverslipRemoved'];
        end
        outputFilePath = fullfile(pathstr, [newFileName '.tiff']);
        yOCT2Tif(modified_data, outputFilePath, 'metadata', metadata, 'clim', clim);
    end

    return;
end

function surface_column = detect_surface_depth(slice, z_size)
    x_size = size(slice, 2);
    % Conditional starting depth based on z_size
    if z_size > 800
        start_depth = 170;
    else
        start_depth = 1;
    end
    
    surface_column = NaN(x_size, 1);  % Initialize surface column with NaNs

    % Detect surface and initially assign NaN for non-confirmable surfaces
    for x = 1:x_size
        for z = start_depth:z_size  % Start from calculated start_depth
            if slice(z, x) >= 2
                if z + 9 <= z_size && all(slice(z+1:z+9, x) > 0)
                    surface_column(x) = z;  % Confirm surface
                    break;
                end
            end
        end
    end

    % Fill in NaN values by interpolating nearest confirmed depths
    for x = 1:x_size
        if isnan(surface_column(x))
            % Find the nearest non-NaN surface depth
            valid_indices = find(~isnan(surface_column));
            if ~isempty(valid_indices)
                % Find the closest valid index
                [~, idx] = min(abs(valid_indices - x));
                closest_valid_index = valid_indices(idx);
                surface_column(x) = surface_column(closest_valid_index);
            end
        end
    end
end

function smoothed_surface = smooth_surface(surface_depth, sigma)
    % Gaussian smoothing across each slice
    kernel_size = ceil(sigma * 3) * 2 + 1;  % Determine the size of the kernel based on sigma
    gaussian_kernel = fspecial('gaussian', [kernel_size, kernel_size], sigma);

    % Replace NaNs with zeros (temporary solution)
    temp_surface = surface_depth;
    nan_mask = isnan(temp_surface);
    temp_surface(nan_mask) = 0;

    % Convolve with the Gaussian kernel, treating zeros as NaNs (this is a simplification)
    smoothed_data = conv2(temp_surface, gaussian_kernel, 'same');

    % Normalize by the effective area of the kernel (excluding NaN contributions)
    normalizing_kernel = conv2(~nan_mask, gaussian_kernel, 'same');
    smoothed_surface = smoothed_data ./ normalizing_kernel;

    % Restore NaNs where the normalizing kernel is zero (no valid data)
    smoothed_surface(normalizing_kernel == 0) = NaN;
    return;
end

function visualization_tiff = create_visualization_tiff(surface_depth, z_size)
    [x_size, y_size] = size(surface_depth);
    max_depth = 500;  % Maximum depth to include in the output
    visualization_tiff = zeros(max_depth, x_size, y_size, 'uint16');  % Using uint16 for higher intensity values
    start_intensity = 40000;  % Start intensity at the surface
    decrease_step = 2000;  % Amount to decrease intensity per depth level

    for x = 1:x_size
        for y = 1:y_size
            z_start = surface_depth(x, y);
            if ~isnan(z_start)
                for z = z_start:min(z_start + start_intensity / decrease_step - 1, max_depth)
                    visualization_tiff(z, x, y) = max(start_intensity - (z - z_start) * decrease_step, 0);
                end
            end
        end
    end
    return;
end

function save_tiff_3d(image_data, filename)
    % Open a new Tiff file for writing
    t = Tiff(filename, 'w');

    % Configure the TIFF tags applicable to all frames
    tagstruct.ImageLength = size(image_data, 1);
    tagstruct.ImageWidth = size(image_data, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16;  % Using 16 bits for higher dynamic range
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = size(image_data, 1);
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression = Tiff.Compression.None;

    % Loop through each slice of image_data
    numSlices = size(image_data, 3);
    for z = 1:numSlices
        t.setTag(tagstruct);
        t.write(image_data(:,:,z));
        if z ~= numSlices
            t.writeDirectory();  % Create a new directory for the next slice
        end
    end

    % Close the Tiff object
    t.close();
end

function modified_data = modify_data_based_on_surface(data, surface_depth)
    % Modify the data based on the detected surface depth
    [z_size, x_size, y_size] = size(data);
    modified_data = data;  % Create a copy to modify

    for x = 1:x_size
        for y = 1:y_size
            surface_z = surface_depth(x, y);
            if ~isnan(surface_z) && surface_z > 1  % Ensure valid surface and avoid edge cases
                % Set all pixels above the surface to NaN
                modified_data(1:surface_z-1, x, y) = NaN;
            end
        end
    end
    disp('-- Data modified based on surface depth.');
end