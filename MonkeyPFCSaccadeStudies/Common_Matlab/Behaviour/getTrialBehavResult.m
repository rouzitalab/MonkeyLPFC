function trials = getTrialBehavResult(trials, ptbParams, analysisParams)
%ptb.eData.trial = getTrialBehavResult(ptb.eData.trial, ptb.params, analysisParams)
%Recalculates trial performance based on detected saccades.
%Checks each saccade's end point against the expected end point.
%New value will be stored in trials(tt).newOutcomeCode;
%analysisParams.critFlip indicates the name of the flipScreen after which
%we search for saccades.

%TODO: critical flipScreen and saccade search window depends on times.

%% Constants and Params

%ptbParams from M&A file.
%                    fixPointSize: 5 or [0 5]
%                       cueHeight: 30 or [0 30]
%                    targetHeight: 30 or [0 30]
%                  targetDistance: 310 < fmap_j_016_00+16
%                  targetDistance: 260 < sra1_j_016_00+10
%                       SR3radius: [0 310] < sra3_2_j_063  px
%               fixPointfixRadius: 75 or [0 90]  px
%            targetPointfixRadius: 75 or [0 125]  px
%                        penWidth: 5 or [0 5]
%  ----- above not present in early fmap files ----
%                          fixRad: 75  %px
%                          fixRad: 100 < fmap_j_002_00+17
%       unitSubjectScreenDistance: 570  ??
%           subjectScreenDistance: 1000 mm
%         subjectScreenResolution: [1024 768] px
%               subjectScreenSize: [735 551] mm
%     subjectScreenPixelPerDegree: 24.4476


%maxSaccadeLatency = 550+300; minSaccadeLatency = 0; %ms
% critFlipScreen = 4; % its actually after 5, but 4 is end of cue

% Get the variables necessary to evaluate saccades based on gaze position.
% trials(x).eyePosition(:, 2:3) are the x,y positions of the
% gaze in units of degrees of visual angle from the top-left of the screen.
% Thus, to easily compare these units to our expectations, we need the
% saccade start point (fixation), saccade end points (targets), and the
% allowable radius for each in units of degrees of visual angle.
% These variables, if given in ptbParams, are in units of pixels.

% Conversion
%Pixel size a.k.a. centimeters per pixel
%pSize = .0702;  % This came from Florian's code.
pSize = norm(ptbParams.subjectScreenSize/10) / norm(ptbParams.subjectScreenResolution);
d2m = ptbParams.subjectScreenDistance/10; %100; %distance to monitor in centimeters
% To convert any pixel measurement to degrees:
% x_in_deg = atand(x_in_px * pSize / d2m);
% x_in_px = d2m * tand(x_in_deg) / pSize;

% Fixation point is in exact center of screen
fixationPx = round(ptbParams.subjectScreenResolution./2);  % Centre of screen, from bottom/top? left, in px
fixationDeg = atand(pSize*fixationPx/d2m);  % Centre of screen, from bottom/top? left, in deg.v.a
%fixationDeg = fixationPx/ptbParams.subjectScreenPixelPerDegree

% Stimuli were squares with length 2*radius
% Get the 'radii' for each stim
%fixation
fixationRadiusDeg = 6;  %4; %deg
if isfield(ptbParams, 'fixPointfixRadius') %&& false
    fixationRadiusPx = ptbParams.fixPointfixRadius(length(ptbParams.fixPointfixRadius));
    fixationRadiusDeg = atand(fixationRadiusPx * pSize / d2m);
end
%target
targetRadiusDeg   = 8; %deg
if isfield(ptbParams, 'targetPointfixRadius')
    targetRadiusPx = ptbParams.targetPointfixRadius(length(ptbParams.targetPointfixRadius));
    targetRadiusDeg = atand(targetRadiusPx * pSize / d2m);
end

clear x y fixationRadiusPx fixationPx locTarg

%% Determine saccade outcomes
% This is better than using the ptb-calculated trial outcomes because
% the saccade data is better, we can relax saccade-timing requirements,
% and we can accommodate some slight shifts in the eyeTracker.

%If we haven't already determined each trial's stimuli, get that now.
if ~isfield(trials, 'targXY')
    trials = getTrialStimInfo(trials, ptbParams);
end

%Set null output.
[trials.newOutcomeCode] = deal(-1);
[trials.sacStartTime] = deal(nan);
[trials.sacEndXY] = deal(nan);
[trials.sacPol] = deal(nan);

%For each trial
nTrials = length(trials);
startPts = nan(nTrials, 2);
endPts = nan(nTrials, 2);
sacTime = nan(nTrials, 1);
sacFlipLims = [find(strcmpi(ptbParams.flipNames, 'targetOnset')) find(strcmpi(ptbParams.flipNames, analysisParams.critFlip))];
for tr_ix = 1:nTrials
    %Get only the first saccade after critical flipScreen
    %TODO: The critical flipScreen may depend on trial type and analysis
    %goals.
    trials(tr_ix).saccades = getSaccadeForTrial(trials(tr_ix), sacFlipLims);
    trials(tr_ix).nSaccades = length(trials(tr_ix).saccades);
    if trials(tr_ix).nSaccades == 1
        startPts(tr_ix, :) = trials(tr_ix).saccades.startPt;
        endPts(tr_ix, :) = trials(tr_ix).saccades.endPt;
        sacTime(tr_ix) = trials(tr_ix).eyePosition(trials(tr_ix).saccades.gazeSampleStart, 1);
    end
end

% Save saccade time to trial structure. Time is in seconds since trial
% start.
sacTime = num2cell(sacTime);
[trials.sacStartTime] = sacTime{:};

%Evaluate fixation OK.
fix_ok = pointInRectangle(startPts, fixationDeg, fixationRadiusDeg);

%Evaluate end point in target, distractor, or neither.
targs = atand(pSize*cat(1, trials.targXY)/d2m);
dists = atand(pSize*cat(1, trials.distXY)/d2m);
targ_dist_bool = false(nTrials, 2);
for tr_ix = 1:nTrials
    if fix_ok(tr_ix)
        targ_dist_bool(tr_ix, :) = pointInRectangle(endPts(tr_ix, :), [targs(tr_ix, :); dists(tr_ix, :)], 2*targetRadiusDeg);
    end
end

% newOutcomeCode is -1. Add 1 (==0) for good saccades, 10 (==9) for
% saccades to distractors.
outcome_add = double(fix_ok & targ_dist_bool(:,1)) + 10*double(fix_ok & targ_dist_bool(:, 2));
new_outcome = num2cell([trials.newOutcomeCode]' + outcome_add);
[trials.newOutcomeCode] = new_outcome{:};

% sacClass is targClass if newOutcomeCode == 0, or distClass if
% newOutcomeCode == 9
targClass = [trials.targClass];
distClass = [trials.targClass];
sacClass = nan(length(trials), 1);
sacClass([trials.newOutcomeCode] == 0) = targClass([trials.newOutcomeCode] == 0);
sacClass([trials.newOutcomeCode] == 9) = distClass([trials.newOutcomeCode] == 9);
sacClass = num2cell(sacClass);
[trials.sacClass] = sacClass{:};

% Get the saccade vector in polar coordinates.
startPtsPx = (d2m/pSize)*tand(startPts);
endPtsPx = (d2m/pSize)*tand(endPts);
sacVecPxXY = endPtsPx - startPtsPx;
[sacTh, sacRad] = cart2pol(sacVecPxXY(:, 1), sacVecPxXY(:, 2)); %sacRad is in cm
sacPol = num2cell([sacTh sacRad], 2);
[trials.sacPol] = sacPol{:};

% Get the saccade endPoint in XY (pixels)
temp = num2cell(endPtsPx, 2);
[trials.sacEndXY] = temp{:};

function boolResult = pointInCircle(pt, circlePt, circleRad)
%function boolResult = pointInCircle(pt, circlePt, circleRad)
nPts = size(pt,1);
nTargs = size(circlePt, 1);
boolResult = false(nPts, nTargs);
for tt = 1:nTargs
    boolResult(:, tt) = sqrt(sum((pt-ones(nPts,1)*circlePt(tt,:)).^2, 2)) <= circleRad;
end

function boolResult = pointInRectangle(pt, rectCent, rectEdgeL)
%function boolResult = pointInRectangle(pt, rectCent, rectEdgeL)
nPts = size(pt,1);
nTargs = size(rectCent, 1);
if length(rectEdgeL) < nTargs
    rectEdgeL = [rectEdgeL repmat(rectEdgeL(1), 1, nTargs-length(rectEdgeL))];
end
boolResult = false(nPts, nTargs);
for tt = 1:nTargs
    rectEdgesLRDU = [...
        rectCent(tt, 1) - rectEdgeL(tt)/2,...
        rectCent(tt, 1) + rectEdgeL(tt)/2,...
        rectCent(tt, 2) - rectEdgeL(tt)/2,...
        rectCent(tt, 2) + rectEdgeL(tt)/2];
    boolResult(:, tt) = ...
        pt(:, 1) >= rectEdgesLRDU(1) &...
        pt(:, 1) <= rectEdgesLRDU(2) &...
        pt(:, 2) >= rectEdgesLRDU(3) &...
        pt(:, 2) <= rectEdgesLRDU(4);
end

function saccade = getSaccadeForTrial(trial, varargin)
global DEF
if nargin < 1
    noSacFlip = find(strcmpi(DEF.flipEvent, 'targetOnset'));  % flipscreen after which no saccade should be allowed until critFlip.
    critFlip = 1 + find(strcmpi(DEF.flipEvent, 'targetOnset'));  % flipscreen after which to search for good saccades.
    if any(strcmpi(trial.newType, 'SR3'))
        critFlip = find(strcmpi(DEF.flipEvent, 'cueOffset'));  % 'fixation offset' more correct, but I'm OK with cue off trials.
    elseif isfield(trial, 'fixPointExtinguishedTime') && ~isempty(trial.fixPointExtinguishedTime)
        critFlip = find(trial.flipScreen == trial.fixPointExtinguishedTime);
    end
else
    noSacFlip = varargin{1}(1);
    critFlip = varargin{1}(2);
end
if length(trial.flipScreen) >= critFlip
    badSamp = find(trial.eyePosition(:,1) >= trial.flipScreen(noSacFlip), 1, 'first');
    critSamp = find(trial.eyePosition(:,1) >= trial.flipScreen(critFlip), 1, 'first');
    sacStarts = [trial.saccades.gazeSampleStart];
    if ~isempty(sacStarts) && ~isempty(critSamp)...
            && any(sacStarts > critSamp) && ~any(sacStarts >= badSamp & sacStarts <= critSamp)
        sac_ix = find(sacStarts > critSamp, 1, 'first');
        saccade = trial.saccades(sac_ix);
    else
        saccade = [];
    end
else
    saccade = [];
end





% % % %% Using correct trials only, determine class->targets (x,y) map
% % % 
% % % expTypes = unique({trials.expType});
% % % expInfo = struct(...
% % %     'type', expTypes,...
% % %     'classes', [],...
% % %     'start_pt', [],...
% % %     'end_pt', [],...
% % %     'vec_theta', [],...
% % %     'vec_amp', []);
% % % 
% % % for exp_ix = 1:length(expTypes)
% % %     
% % %     % Identify correct trials for this experiment type.
% % %     tr_bool = strcmpi({trials.expType}, expTypes{exp_ix});
% % %     tr_bool = tr_bool & [trials.outcomeCode]==0;
% % %     
% % %     %Pre-allocate
% % %     n_trials = sum(tr_bool);
% % %     all_trials = trials(tr_bool);
% % %     start_stop = nan(n_trials, 4);  % x_start y_start x_stop y_stop
% % %     vecTheta = nan(n_trials, 1);
% % %     vecAmp = nan(n_trials, 1);
% % %     % For each good saccade
% % %     for tr_ix = 1:n_trials
% % %         saccade = getSaccadeForTrial(all_trials(tr_ix));  % TODO: Use some time limits.
% % %         start_stop(tr_ix, :) = [saccade.startPt saccade.endPt];
% % %         vecTheta(tr_ix) = saccade.vectDir;
% % %         vecAmp(tr_ix) = saccade.vectAmp;
% % %     end
% % %     % Determine the nearest direction of possible directions (targetTheta)
% % %     vecTheta(vecTheta<-45/2) = vecTheta(vecTheta<-45/2)+360;
% % %     
% % %     % Start point is class independent.
% % %     expInfo(exp_ix).start_pt = round(10*mean(start_stop(:, 1:2)))./10;
% % %     expInfo(exp_ix).start_pt_exp = fixationDeg;
% % %     
% % %     % Fill in the per-class information.
% % %     expInfo(exp_ix).classes = unique([trials(tr_bool).class]);
% % %     n_classes = length(expInfo(exp_ix).classes);
% % %     expInfo(exp_ix).end_pt = nan(n_classes, 2);
% % %     expInfo(exp_ix).end_pt_exp = nan(n_classes, 2);
% % %     expInfo(exp_ix).end_pt_dif = nan(n_classes, 1);
% % %     expInfo(exp_ix).vec_theta = nan(n_classes, 1);
% % %     expInfo(exp_ix).vec_amp = nan(n_classes, 1);
% % %     for cl_ix = 1:n_classes
% % %         
% % %         % Identify trials in this class.
% % %         cl_bool = [all_trials.class] == expInfo(exp_ix).classes(cl_ix);
% % %         
% % %         % End point
% % %         expInfo(exp_ix).end_pt(cl_ix, :) = round(10*mean(start_stop(cl_bool, 3:4)))/10;  % average
% % %         %nearest target
% % %         targ_diff = sqrt(sum((targetPtsDeg - ones(size(targetPtsDeg,1),1)*expInfo(exp_ix).end_pt(cl_ix,:)).^2, 2));
% % %         match_bool = targ_diff == min(targ_diff);
% % %         expInfo(exp_ix).end_pt_exp(cl_ix, :) = targetPtsDeg(match_bool, :);
% % %         expInfo(exp_ix).end_pt_dif(cl_ix) = targ_diff(match_bool);
% % %         
% % %         % Vector polar measurements.
% % %         expInfo(exp_ix).vec_theta(cl_ix) = mean(vecTheta(cl_bool));
% % %         expInfo(exp_ix).vec_amp(cl_ix) = mean(vecAmp(cl_bool));
% % %     end
% % % end
% % % ptb.eData.expInfo = expInfo;
% % % clear exp_ix expTypes tr_bool n_trials all_trials start_stop tr_ix saccade
% % % clear  n_classes cl_ix cl_bool targ_diff match_bool