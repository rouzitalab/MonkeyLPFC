function gazeData = invalidateGazeData(gazeData)
%function gazeData = invalidateGazeData(gazeData)
%gazeData.state = whether each sample is valid or not.
global DEF

%Define a limit rectangle 40deg lr and 30 deg ud +10%
visualAreaDegLimit.up   = 30+3; %degrees from lowerLeft
visualAreaDegLimit.down = 0-3;
visualAreaDegLimit.left = 0-4;
visualAreaDegLimit.right= 40+4;
lowPupilDelay = 5; %samples. Buffer to invalidate Datapoints before and after blink
minValidBlockDur = 500; %samples. Only consider sample good if in a block at least this long.

gazeData.state = zeros(size(gazeData.pupil),DEF.errorcodes)+DEF.EYE_State_undefined;

% Offscreen
offBool = gazeData.degree(:,1) > visualAreaDegLimit.right | ...
    gazeData.degree(:,1) < visualAreaDegLimit.left | ...
    gazeData.degree(:,2) > visualAreaDegLimit.up | ...
    gazeData.degree(:,2) < visualAreaDegLimit.down;
gazeData.state(offBool) = DEF.EYE_State_offscreen;

% Blinks. Overwrites offscreen state.
closedBool = gazeData.pupil==0;
blinkEdges = find(diff([closedBool(1); closedBool])~=0);
%buffer around blinkEdges
for be_ix = 1:length(blinkEdges)
    start_ix = max([1 blinkEdges(be_ix)-lowPupilDelay]);
    stop_ix = min([blinkEdges(be_ix)+lowPupilDelay length(closedBool)]);
    closedBool(start_ix : stop_ix) = true;
end
gazeData.state(closedBool) = DEF.EYE_State_invalid;

% nans. The asc data files actually have nans.
gazeData.state(isnan(gazeData.gaze)) = DEF.EYE_State_invalid;

% Only include blocks of valid samples at least minValidBlockDur long.
yetUndefBool = gazeData.state == DEF.EYE_State_undefined;
yetUndefStarts = find(diff(yetUndefBool)==1);
if yetUndefBool(1)
    yetUndefStarts = [true;yetUndefStarts];
end
yetUndefStops = find(diff(yetUndefBool)==-1);
if yetUndefBool(end)
    yetUndefStops = [yetUndefStops;length(yetUndefBool)];
end   

badBlocksBool = (yetUndefStops-yetUndefStarts) < minValidBlockDur;
badBlockStarts = yetUndefStarts(badBlocksBool);
badBlockStops = yetUndefStops(badBlocksBool);

for block_ix = 1:length(badBlockStarts)
    gazeData.state(badBlockStarts(block_ix):badBlockStops(block_ix)) = DEF.EYE_State_invalid;
end