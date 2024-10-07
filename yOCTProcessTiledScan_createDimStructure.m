function [dimOneTile, dimOutput] = yOCTProcessTiledScan_createDimStructure(tiledScanInputFolder)
% This is an auxilary function of yOCTProcessTiledScan designed to retun a
% dim structure for both a single tile and the entier tiledScan

%% Search and Load JSON file from  the tiledScanInputFolder
json = awsReadJSON([tiledScanInputFolder 'ScanInfo.json']);

%% Pretend processing the first scan to get the dimension structure of one tile
firstDataFolder = [tiledScanInputFolder json.octFolders{1}];
if (not(awsExist(firstDataFolder,'dir')))
    error('%s does not exist.', firstDataFolder);
end

dimOneTile = yOCTLoadInterfFromFile(firstDataFolder,'OCTSystem',json.OCTSystem,'peakOnly',true);
tmp = zeros(size(dimOneTile.lambda.values(:)));
dimOneTile = yOCTInterfToScanCpx (tmp, dimOneTile, 'n', json.tissueRefractiveIndex, 'peakOnly',true);

% Update dimensions to mm
dimOneTile = yOCTChangeDimensionsStructureUnits(dimOneTile, 'mm');

%% Update X&Y positions as they might not be reliable from scan

% Dimensions of one tile
if ~isfield(json,'xRange_mm')
    % Backward compatibility
    warning('Note, that "%s" contains an old version scan, this will be depricated by Jan 1st, 2025',tiledScanInputFolder)
    dimOneTile.x.values = json.xOffset+json.xRange*linspace(-0.5,0.5,json.nXPixels+1);
    dimOneTile.y.values = json.yOffset+json.yRange*linspace(-0.5,0.5,json.nYPixels+1);
else
    dimOneTile.x.values = json.xOffset+json.tileRangeX_mm*linspace(-0.5,0.5,json.nXPixels+1);
    dimOneTile.y.values = json.yOffset+json.tileRangeY_mm*linspace(-0.5,0.5,json.nYPixels+1);
end
dimOneTile.x.values(end) = [];
dimOneTile.y.values(end) = [];

%% Compute pixel size
dx = diff(dimOneTile.x.values(1:2));
if (length(dimOneTile.y.values) > 1)
    dy = diff(dimOneTile.y.values(1:2));
else
    dy = 0; % No y axis
end
dz = diff(dimOneTile.z.values(1:2));

%% Compute dimensions of the entire output

zDepths_mm = json.zDepths;
if ~isfield(json,'xCenters_mm')
    % Backward compatibility
    xCenters_mm = json.xCenters;
    yCenters_mm = json.yCenters;
else
    xCenters_mm = json.xCenters_mm;
    yCenters_mm = json.yCenters_mm;
end

xAll_mm = (min(xCenters_mm)+dimOneTile.x.values(1)):dx:(max(xCenters_mm)+dimOneTile.x.values(end)+dx);xAll_mm = xAll_mm(:);xAll_mm = xAll_mm.';
yAll_mm = (min(yCenters_mm)+dimOneTile.y.values(1)):dy:(max(yCenters_mm)+dimOneTile.y.values(end)+dy);yAll_mm = yAll_mm(:);yAll_mm = yAll_mm.';
zAll_mm = (min(zDepths_mm )+dimOneTile.z.values(1)):dz:(max(zDepths_mm) +dimOneTile.z.values(end));zAll_mm = zAll_mm(:);

% Correct for the case of only one scan
if (length(xCenters_mm) == 1) %#ok<ISCL>
    xAll_mm = dimOneTile.x.values;
end
if (length(yCenters_mm) == 1) %#ok<ISCL> 
    yAll_mm = dimOneTile.y.values;
end

%% Create dimensions data structure
dimOutput.lambda = dimOneTile.lambda;
dimOutput.z = dimOneTile.z; % Template, we will update it soon
dimOutput.z.values = zAll_mm(:)' - zAll_mm(1);
dimOutput.z.origin = 'z=0 is the top of OCT image when zDepths=0 scan was taken';
dimOutput.x = dimOneTile.x;
dimOutput.x.origin = 'x=0 is OCT scanner origin when xCenters=0 scan was taken';
dimOutput.x.values = xAll_mm(:)';
dimOutput.x.index = 1:length(dimOutput.x.values);
dimOutput.y = dimOneTile.y;
dimOutput.y.values = yAll_mm(:)';
dimOutput.y.index = 1:length(dimOutput.y.values);
dimOutput.y.origin = 'y=0 is OCT scanner origin when yCenters=0 scan was taken';
dimOutput.aux = dimOneTile.aux;


