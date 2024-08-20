% Run this script to detect the surface of the tissue from a file, generate
% a heatmap and create a file with a surface representation

function process_oct_data()
    % Load the full z-stack volume
    filePath = 'C:\Users\Alejo\OneDrive\Escritorio\Old Code\OCTHist Code - Git Old Working File\myOCT\Diego_Files\August_Tests\3D_Volume_Reconstructions\1um_Large_Scan\Data08.tif';
    [data, metadata, c] = yOCTFromTif(filePath);
    disp([newline, '-- Data loaded successfully.']);

    % Extract directory, base name, and extension for the input file
    [fileDir, baseName, ~] = fileparts(filePath);
    outputBaseName = [baseName '_surface'];
    outputFileExt = '.tif';
    outputFileName = fullfile(fileDir, [outputBaseName outputFileExt]);
    counter = 0;

    % Check if file exists, and create a new file name if necessary
    while exist(outputFileName, 'file')
        counter = counter + 1;
        outputFileName = fullfile(fileDir, [outputBaseName num2str(counter, '%02d') outputFileExt]);
    end

    % Dimensions of the data
    [z_size, x_size, y_size] = size(data);

    % Initialize an array to hold the surface depth for each (x, y) coordinate
    surface_depth = zeros(x_size, y_size);

    % Process each slice
    for y = 1:y_size
        % Extract the current slice
        current_slice = data(:, :, y);

        % Process the slice to find the surface depth
        surface_depth(:, y) = find_surface(current_slice, z_size);
    end

    % Apply Gaussian smoothing to the surface depth data
    sigma = 1.5;  % Standard deviation for the Gaussian kernel
    smoothed_surface_depth = gaussian_smooth_surface(surface_depth, sigma);

    % Round the smoothed surface depth to the nearest integer
    smoothed_surface_depth = round(smoothed_surface_depth);

    % Display the smoothed and rounded surface depth as a heatmap
    figure;
    imagesc(smoothed_surface_depth.');
    colormap(flipud(jet));
    colorbar;
    set(gca, 'YDir', 'normal');
    title('Surface Depth Heatmap');
    xlabel('X-axis');
    ylabel('Y-axis');

    % Create a 3D visualization TIFF
    visualization_tiff = create_visualization_tiff(smoothed_surface_depth, z_size);

    % Save the 3D visualization
    save_tiff_3d(visualization_tiff, outputFileName);

    disp([newline, '-- Visualization TIFF created and saved at ', outputFileName, ' --', newline, '-- Process completed.']);
end

function surface_column = find_surface(slice, z_size)
    x_size = size(slice, 2);
    % Conditional starting depth based on z_size
    if z_size > 800
        start_depth = 60;
    else
        start_depth = 1;
    end
    
    surface_column = NaN(x_size, 1);  % Initialize surface column with NaNs

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
    return;
end

function smoothed_surface = gaussian_smooth_surface(surface_depth, sigma)
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
    tagstruct.Software = 'MATLAB';
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
    return;
end
