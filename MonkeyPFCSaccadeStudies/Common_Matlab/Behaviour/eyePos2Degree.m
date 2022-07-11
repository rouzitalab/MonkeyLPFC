function gazeData = eyePos2Degree(gazeData, lineOfSightPix, minPixPerDegree, minViewDist)
%function gazeData = eyePos2Degree(gazeData, lineOfSightPix, minPixPerDegree, minViewDist)
%I have not looked this over.
global DEF CONST
gazeData.degree = nan(size(gazeData.gaze), DEF.eyePosDegreeFormat);
if ~isfield(gazeData,'resolution')
    gazeData.resolution = nan(size(gazeData.gaze), DEF.eyePosResolutionFormat);
end

minCentimeterPerDegree = minViewDist/CONST.dist1cmPerDeg;
pixelSize              = 1/(minPixPerDegree/minCentimeterPerDegree);
minViewPixelDist       = minViewDist/pixelSize;
upperLeft2CenterDist(1)= atand(lineOfSightPix(1)/minViewPixelDist);
upperLeft2CenterDist(2)= atand(lineOfSightPix(2)/minViewPixelDist);

% Resolution
gazeData.resolution(:,1) = sqrt(minViewDist^2 + (pixelSize.*(gazeData.gaze(:,1)-lineOfSightPix(1))).^2) / CONST.dist1cmPerDeg/pixelSize;
gazeData.resolution(:,2) = sqrt(minViewDist^2 + (pixelSize.*(gazeData.gaze(:,2)-lineOfSightPix(2))).^2) / CONST.dist1cmPerDeg/pixelSize;

% degree: corrected for centerPixel
% Inverse tangent; result in degrees
gazeData.degree(:,1) = atand( (gazeData.gaze(:,1) - lineOfSightPix(1)) /minViewPixelDist) + upperLeft2CenterDist(1);
gazeData.degree(:,2) = atand( (gazeData.gaze(:,2) - lineOfSightPix(2)) /minViewPixelDist) + upperLeft2CenterDist(2);