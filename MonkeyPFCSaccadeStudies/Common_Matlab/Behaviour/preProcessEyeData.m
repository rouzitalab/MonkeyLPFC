function gazeData = preProcessEyeData(gazeData, varargin)
fprintf('Preprocessing gazeData... ');

% Default parameters
defParams.adjustments = [];
defParams.gazeFilter = struct('function', 'runavg', 'param', 3);
defParams.ptbParams = struct('subjectScreenSize', 2*[512 384], ...
    'subjectScreenPixelPerDegree', 25);
defParams.minViewDist = 100;
global DEF

% Merge varargin with defaults.
params = varg2params(varargin, defParams);

% Merging may have wiped out some defaults without replacing them.
if ~isfield(params.ptbParams, 'subjectScreenSize')
    params.ptbParams.subjectScreenSize = defParams.ptbParams.subjectScreenSize;
end
if ~isfield(params.ptbParams, 'subjectScreenPixelPerDegree')
    params.ptbParams.subjectScreenPixelPerDegree = defParams.ptbParams.subjectScreenPixelPerDegree;
end

% correct for online adjustments in eData.eyePosCalibAdjust
gazeData = correctForEyePosAdjustments(gazeData, params.adjustments);

% lowpass-filter gazedata.gaze
if ~isempty(varargin)
    gazeFilterKernel = makeGazeFilterKernel(params.gazeFilter);
    gazeData.gaze(:,1) = conv2(gazeData.gaze(:,1), gazeFilterKernel', 'same');
    gazeData.gaze(:,2) = conv2(gazeData.gaze(:,2), gazeFilterKernel', 'same');
end

% convert the PixelPos in degree & pixel resolution
lineOfSightPix = params.ptbParams.subjectScreenSize/2;
minPixPerDegree = params.ptbParams.subjectScreenPixelPerDegree;
gazeData = eyePos2Degree(gazeData, lineOfSightPix, minPixPerDegree, params.minViewDist);

% define limits of degrees beyond which we assume the eyeposition to be 
gazeData = invalidateGazeData(gazeData);

% Calculate velocity, using sample times.
gazeData.velocity               = nan(size(gazeData.gaze));
gazeData.velocity(2:end,:)      = bsxfun(@rdivide, diff(double(gazeData.degree)), diff(double(gazeData.sample))); % degree displacement between samples / time step between samples
gazeData.velocity               = cast(gazeData.velocity, DEF.eyePosVelocityFormat);
gazeData.velUnits               = 'degrees/msec';
gazeData.acceleration           = nan(size(gazeData.gaze), DEF.eyePosVelocityFormat);
gazeData.acceleration(2:end,:)  = diff(gazeData.velocity); % degree displacement between samples
gazeData.accelUnits             = 'degrees/msec^2';
gazeData.dir                    = rad2deg(atan2(gazeData.velocity(:,2), gazeData.velocity(:,1))); 
gazeData.velNorm                = sqrt(gazeData.velocity(:,1).^2 + gazeData.velocity(:,2).^2);
%In Florian's code, speed was sqrt( (eyeData.sampleRate*vel1).^2 + (eyeData.sampleRate*vel2).^2 )
%In other words, vel is first converted to deg/sec.
%gazeData.speed                  = sqrt(gazeData.velocity(:,1).^2 + gazeData.velocity(:,2).^2);
%This creates an inconsistency in fields' units so I got rid of it.
fprintf('Done.\n');