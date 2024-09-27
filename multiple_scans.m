% This is the same as Demo_ScanAndProcess_3D Version 9.27.2024
% BUT this will perform multiple scans according to input

% Please, start with the Multiple Scans inputs, and continue with the
% other inputs as usual.

%% Inputs
% Multiple Scans Inputs

% Set the number of scans and the interval between them
numScans = 5; % Number of scans
intervalMinutes = 10; % Interval in minutes between scans

% Define the 3D Volume and scanning parameters

pixel_size_um = 1; % x-y Pixel size in microns
xOverall_mm = [-0.25 0.25]; % Overall volume to scan [start, finish]
yOverall_mm = [-0.25 0.25]; % Overall volume to scan [start, finish]

octProbePath = yOCTGetProbeIniPath('40x','OCTP900'); % Probe ini path
octProbeFOV_mm = 0.5; % Field of view from the probe
oct2stageXYAngleDeg = 0; % Angle between x axis of the motor and the Galvo's x axis

scanZJump_um = 5; % Z stack increment
zToScan_mm = ((-40:scanZJump_um:250))*1e-3; %[mm]
focusSigma = 20; % Size of each focus in z [pixel]

tissueRefractiveIndex = 1.33; % Refractive index of the tissue
dispersionQuadraticTerm = -1.466e+08; % Dispersion term

output_folder = '\'; % Directory to save scan files

skipScanning = false; % Set to true to skip scanning
skipProcessing = false; % Set to true to skip processing

focusPositionInImageZpix = [393]; % Known depth of focus position

%% Execute multiple scans

for scanIndex = 1:numScans

    currentStartTime = datetime('now');
    fprintf('\n\nStarting scan %d/%d at %s...\n\n', scanIndex, numScans, datestr(currentStartTime, 'HH:MM AM'));
    
    perform_scan_and_process(scanIndex, output_folder, pixel_size_um, xOverall_mm, yOverall_mm, octProbePath, octProbeFOV_mm, oct2stageXYAngleDeg, zToScan_mm, focusSigma, tissueRefractiveIndex, dispersionQuadraticTerm, skipScanning, skipProcessing, focusPositionInImageZpix);
    
    currentEndTime = datetime('now');
    fprintf('\n\nScan %d/%d finished at %s...\n', scanIndex, numScans, datestr(currentEndTime, 'HH:MM AM'));

    if scanIndex < numScans
        fprintf('\nNext scan will start in %d minutes at %s.\n', intervalMinutes, datestr(nextStartTime, 'HH:MM AM'));
        fprintf('\nSafe to pause now. \nOnce the next scan starts, stopping is not advised.\n\n');
        pause(intervalMinutes * 60); % Pause between scans in secs.
    end
    
end

fprintf('\nAll requested scans were completed.\n');

function perform_scan_and_process(scanIndex, output_folder, pixel_size_um, xOverall_mm, yOverall_mm, octProbePath, octProbeFOV_mm, oct2stageXYAngleDeg, zToScan_mm, focusSigma, tissueRefractiveIndex, dispersionQuadraticTerm, skipScanning, skipProcessing, focusPositionInImageZpix)
    %% Compute scanning parameters

    % Check that sufficient ammount of gel is above the tissue for proper focus
    if (min(zToScan_mm)) > -100e-3
        warning('Because we use gel above tissue to find focus position. It is important to have at least one of the z-stacks in the gel. Consider having the minimum zToScan_mm to be -100e-3[mm]')
    end

    %% Create unique folder for each scan
    scanFolder = sprintf('%s/Scan_%d', output_folder, scanIndex);

    if ~exist(scanFolder, 'dir')
        mkdir(scanFolder);  % Create the folder if it doesn't exist
    end
    volumeOutputFolder = [scanFolder '/OCTVolume'];

    if ~exist(volumeOutputFolder, 'dir')
        mkdir(volumeOutputFolder);  % Assures the volume directory is created
    end

    %% Perform the scan
    disp('Please adjust sample such that the sample-gel interface is at OCT focus')
    
    fprintf('%s Scanning Volume in %s\n', datestr(datetime), volumeOutputFolder);
    scanParameters = yOCTScanTile (...
        volumeOutputFolder, ...
        xOverall_mm, ...
        yOverall_mm, ...
        'octProbePath', octProbePath, ...
        'tissueRefractiveIndex', tissueRefractiveIndex, ...
        'octProbeFOV_mm', octProbeFOV_mm, ...
        'pixelSize_um', pixel_size_um, ...
        'xOffset',   0, ...
        'yOffset',   0, ... 
        'zDepths',   zToScan_mm, ... [mm]
        'oct2stageXYAngleDeg', oct2stageXYAngleDeg, ...
        'skipHardware',skipScanning, ...
        'v',true  ...
        );

    %% Find focus in the scan
    if isempty(focusPositionInImageZpix)
        fprintf('%s Find focus position volume\n',datestr(datetime));
        focusPositionInImageZpix = yOCTFindFocusTilledScan(volumeOutputFolder,...
            'reconstructConfig',{'dispersionQuadraticTerm',dispersionQuadraticTerm},'verbose',true);
    end

    %% Process the scan
    if ~skipProcessing
        if exist(volumeOutputFolder, 'dir')  % Check if the scan data exists
            enableCropping = true;  % Set to true to enable cropping, false to disable
            fprintf('%s Processing\n', datestr(datetime));
            outputTiffFile = sprintf('%s/Recons_Scan_%d.tiff', scanFolder, scanIndex);
            yOCTProcessTiledScan(...
                volumeOutputFolder, ... Input
                {outputTiffFile},... Save only Tiff file as folder will be generated after smoothing
                'focusPositionInImageZpix', focusPositionInImageZpix,... No Z scan filtering
                'focusSigma',focusSigma,...
                'dispersionQuadraticTerm',dispersionQuadraticTerm,... Use default
                'cropZAroundFocusArea', enableCropping, ...
                'interpMethod','sinc5', ...
                'v',true);
        else
            fprintf('%s Processing skipped as no scan data available\n', datestr(datetime));
        end
    else
        fprintf('%s Processing is disabled. Change skipProcessing to false to enable.\n', datestr(datetime));
    end
end
