% Run this demo to use Thorlabs system to scan a 3D OCT Volume and process
% it.
% Before running this script, make sure myOCT folder is in path for example
% by running: addpath(genpath('F:\Jenkins\Scan OCTHist Dev\workspace\'))

% The protocol for how to use this script can be found here:
% https://docs.google.com/presentation/d/1EUYneJwzGAgj2Qg-rG0k6EQb5t1KOCVi0VxEJl8mPmM/edit#slide=id.g25bcdbd2c45_0_0

%% Inputs

% Define the 3D Volume
pixel_size_um = 1; % x-y Pixel size in microns
xOverall_mm = [-0.25 0.25]; % Define the overall volume you would like to scan [start, finish]. For 10x use [-0.5 0.5] for 40x use [-0.25 0.25]
yOverall_mm = [-0.25 0.25]; % Define the overall volume you would like to scan [start, finish]. For 10x use [-0.5 0.5] for 40x use [-0.25 0.25]
% Uncomment below to scan one B-Scan.
% yOverall_mm = 0;

% Define probe 
octProbePath = yOCTGetProbeIniPath('40x','OCTP900'); % Probe ini spec, you can use yOCTGetProbeIniPath('10x','OCTP900') etc
octProbeFOV_mm = 0.5; % How much of the field of view to use from the probe.
oct2stageXYAngleDeg = 0; % Angle between x axis of the motor and the Galvo's x axis

% Define z stack and z-stitching
scanZJump_um = 5; % Use 15 microns for 10x lens, 5 microns for 40x lens
zToScan_mm = ((-40:scanZJump_um:250))*1e-3; %[mm]
focusSigma = 20; %When stitching along Z axis (multiple focus points), what is the size of each focus in z [pixel], use 20 for 10x, 1 for 40x

% Other scanning parameters
tissueRefractiveIndex = 1.33; % Use either 1.33 or 1.4 depending on the results. Use 1.4 for brain.
%dispersionQuadraticTerm=6.539e07; % 10x
%dispersionQuadraticTerm=9.56e7;   % 40x
dispersionQuadraticTerm=-1.466e+08;  % 10x, OCTP900

% Where to save scan files
output_folder = '\';

% Set to true if you would like to skip either Scanning or Processing
skipScanning = true;     % Skip scanning if set to true
skipProcessing = false;  % Skip processing if set to true

% If depth of focus position is known, write it here. If you would like the script to help you keep empty
focusPositionInImageZpix = [393];

%% Compute scanning parameters

% Check that sufficient ammount of gel is above the tissue for proper focus
if (min(zToScan_mm)) > -100e-3
    warning('Because we use gel above tissue to find focus position. It is important to have at least one of the z-stacks in the gel. Consider having the minimum zToScan_mm to be -100e-3[mm]')
end

%% Perform the scan
volumeOutputFolder = [output_folder '/OCTVolume/'];
disp('Please adjust sample such that the sample-gel interface is at OCT focus')

fprintf('%s Scanning Volume\n',datestr(datetime));
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
        enableCropping = false;  % Set to true to enable cropping, false to disable
        fprintf('%s Processing\n', datestr(datetime));
        outputTiffFile = [output_folder '/Image.tiff'];
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
