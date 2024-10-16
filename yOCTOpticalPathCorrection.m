function yOCTOpticalPathCorrection(opticalPathCorrectionOptions, outputPath)
% This function performs optical path correction on a 3D TIFF file and saves the corrected image.
%
% INPUTS:
%   - outputTiffFile: Path to the TIFF file containing the 3D image data with dimensions (z, x, y).
%   - opticalPathCorrectionOptions: This variable can be one of the following options:
%       1. Path to the Volumes folder where the ScanInfo.json file is stored.
%       2. Struct containing information extracted from ScanInfo.json.
%       3. Array of polynomial coefficients used for optical path correction, must contain 12 elements.
%   - outputPath: Path where the final corrected TIFF file will be saved.

    % Extract optical path polynomial coefficients from opticalPathCorrectionOptions
    if isstruct(opticalPathCorrectionOptions)
        OP_p = opticalPathCorrectionOptions.octProbe.OpticalPathCorrectionPolynomial;
        OP_p = OP_p(:)';
    elseif ischar(opticalPathCorrectionOptions) || isstring(opticalPathCorrectionOptions)
        inputVolumeFolder = awsModifyPathForCompetability([opticalPathCorrectionOptions '/']);
        json = awsReadJSON([inputVolumeFolder 'ScanInfo.json']);
        OP_p = json.octProbe.OpticalPathCorrectionPolynomial;
        OP_p = OP_p(:)';
    elseif isnumeric(opticalPathCorrectionOptions)
        if ~(isequal(size(opticalPathCorrectionOptions), [1,12]) || isequal(size(opticalPathCorrectionOptions), [12,1]))
            error('Optical path correction polynomial must have 12 terms')
        end
        OP_p = opticalPathCorrectionOptions(:)';
    else
        error('opticalPathCorrectionOptions must be a file path to the ScanInfo.json file, a struct representing the json, or an array of polynomial terms.');
    end
    
    % Fix path
    outputPath = char(awsModifyPathForCompetability(outputPath));
  
    % Separate the coefficients for x and y
    pX = OP_p(1:6);
    pY = OP_p(7:12);
    
    % Load TIF data
    [data, ~, ~] = yOCTFromTif(outputPath);
    [zSize, xSize, ySize] = size(data);
    
    % Prepare data based on polynomios
    [X, Y] = meshgrid(1:xSize, 1:ySize);
    Z_pred = polyval(pX, Y) + polyval(pY, X); % Combine polynomios
    
    % Apply correction
    correctedData = zeros(size(data), 'like', data);
    
    for x = 1:xSize
        for y = 1:ySize
            zShift = round(Z_pred(y, x) - Z_pred(250, 250)); % Based on centers
            correctedData(:, x, y) = circshift(data(:, x, y), -zShift); % Adjust using circshift
        end
    end
    
    % Save file
    yOCT2Tif(correctedData, outputPath);
    fprintf('Optical path correction has been applied and the corrected file is saved.\n');

end
