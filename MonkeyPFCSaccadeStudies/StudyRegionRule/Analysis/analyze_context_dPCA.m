% Do 4 dPCAs for RuleA vs RuleB (U v L, U v R, D v L, D v R)
% Cue-colour will be other dimension in dPCA when available.
% For now I am ignoring the location of the distractor.

%% Paths and Constants
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths
sessions = my_sessions('RegionRule');
analysisParams = my_anaparams('dPCA');
%addpath(fullfile(paths.mlDevel, 'Helper'));
addpath(genpath(fullfile(paths.ml3rd, 'DataHigh')));
addpath(fullfile(paths.ml3rd, 'dPCA', 'matlab'));

%% Constants and other frequently accessed variables

%params struct for DataHigh.
dhParams = struct(...
        'binWidth', analysisParams.binWidth,...
        'use_sqrt', true,...
        'kern', analysisParams.kernSD,...
        'trial_average', false,...
        'keep_neurons', nan);  % true(nUnits, 1);

nAna = length(analysisParams.anaWins);
nSessions = length(sessions);



%% OLD
% % Relate saccade direction to a class by binning. Class order is arbitrary.
% % Note that the eye-tracker origin is top-left of screen, so y +ve is down.
% classTargCentres = -90:45:225;
% classTargCentres(classTargCentres<0) = classTargCentres(classTargCentres<0) + 360;
% classStr = {'UU' 'UR' 'RR' 'DR' 'DD' 'DL' 'LL' 'UL'};
% nClasses = length(classTargCentres);

%% Init variables we will save.
sess_unit_id = cell(nSessions, 1);  % each cell is nNeurons x 2 = [sess_ix unit_ix]
trial_class = cell(nSessions, 1);  % each cell is nClasses x 1
trialNum = cell(nSessions, 1);
firingRates = cell(nSessions, nAna);  % each cell is nNeurons x nTimes x nClasses x nTrials

%nTimes and nClasses will have to be permuted prior to dPCA

sess_ix = 4;
%%
for sess_ix = 1:nSessions
    %%
    this_sess = sessions(sess_ix);
    
    %% Load and process ptbmat
    ptb = load(fullfile(paths.preprocessed, 'ptb', this_sess.ptbfname), '-mat');
    [ptb.eData.trial, ptb.params] = getTrialStimInfo(ptb.eData.trial, ptb.params);
    ptb.eData.trial = getTrialBehavResult(ptb.eData.trial, ptb.params, analysisParams);
    [IsInCondition] = getCornerConditionMembership(ptb.eData.trial);
    %squeeze(sum(any(IsInCondition, 4)))
    ptb.eData.trial = getNewClass(ptb.eData.trial, analysisParams);
    
    %% Load and process cerbDat
    cDat = load(fullfile(paths.preprocessed, 'cDat', [this_sess.edfname(1:end-3) 'mat']));
    
    % Consolidate cerbDat ptbmat
    cDat = consolidate_ptb_cdat(ptb, cDat);
    
    % Eliminate trials
    % - without saccades
    % - without timeLockEvent
    % - without appropriate newOutcomeCodes
    % - without minTrialsPerGroup
    trialBool = triageTrials(cDat, analysisParams);
    cDat.trial = cDat.trial(trialBool);
    
%     % Further eliminate trials without enough samples.
%     cDat = trimTrials(cDat, analysisParams.anaWins(end));
%     aftSac = round(cellfun(@length, {cDat.trial.tVec}) - [cDat.trial.sacStartTime]);
%     trialBool = trialBool & aftSac >= analysisParams.anaWins(end).winEdges(end);
%     cDat.trial = cDat.trial(trialBool); clear trialBool aftSac
    
    % Eliminate units that fire < 1Hz
    % Adds cDat.invalidUnits.
    % If thisParams.useUnits == 'merged' then it merges all units per
    % channel into one unit per channel.
    cDat = triageUnits(cDat, analysisParams);
    
    %% Easy-access to classes
    classId = [cDat.trial.newClass];
    [uqClasses, inClassIx] = unique(classId);
    classCount = hist(classId, uqClasses);
    nClasses = length(uqClasses);
    nTrials = length(cDat.trial);
    nUnits = size(cDat.trial(1).raster, 2);
    dhParams.keep_neurons = true(nUnits, 1);
    clear inClassIx
    
    %% Prepare data for dPCA
    
    % Identify eventTime, avoidTime, a tVec, and timeBool for each trial.
    cDat = trimTrials(cDat, analysisParams.anaWins(1));
    
    % Convert trial data to data structure expected by DataHigh
    D = cDat2DataHigh(cDat.trial);
    
    % Reduce the dimensionality by binning and smoothing spikes.
    % Change nan to a different algorithm for PCA, PPCA, FA, LDA, GPFA
    [newD, ~, ~] = reducedims(D, -1, nUnits, dhParams);
    
    % Our data lengths may be unequal. Pad with nans.
    % Note that epochStarts will still be broken.
    dSize = nan(nTrials, 1);
    for tr_ix = 1:nTrials
        dSize(tr_ix) = size(newD(tr_ix).data, 2);
    end
    nBins = max(dSize);
    for tr_ix = 1:nTrials
        newD(tr_ix).data = [newD(tr_ix).data nan(nUnits, nBins - dSize(tr_ix))];
    end
    
    %Get the time vector
    tr_ix = find(dSize == max(dSize), 1, 'first');
    tVec = cDat.trial(tr_ix).tVec(cDat.trial(tr_ix).timeBool);
    tVec = tVec(1, 1:dhParams.binWidth * floor(size(tVec, 2)/dhParams.binWidth));  % Cut off excess
    tVec = reshape(tVec, dhParams.binWidth, []);
    tVec = mean(tVec);
    
    %datahigh2dpca(newD)
    newD = cat(3, newD.data);
    nTimes = size(newD, 2);
    
    clear D tr_ix
    
    % For each corner
    for cp = 1:size(corner_pairs, 1)
        
        % Identify the trials that belong to each rule & colour for this
        % corner
        corner_bool = false(nTrials, 2, 3);  % 2 in corner, 3 cue colours
        for ix = 1:2
            for colour_ix = 1:3
                corner_bool(:, ix, colour_ix) = classId == (corner_pairs(cp, ix)+colour_ix-1);
            end
        end
        
        %We will only analyze rule x colour combinations for which there
        %are trials.
        good_colour = find(squeeze(~sum((sum(corner_bool, 1) == 0), 2)));
        max_trials = max(max(squeeze(sum(corner_bool(:, :, good_colour)))));
        nColours = length(good_colour);
        
        %Allocate memory
        %parameters: 1 - context, 2 - colour, 3 - time
        firingRates = nan(nUnits, 2, nColours, nTimes, max_trials);
        
        %Get the trial data into the 5-D matrix
        for colour_ix = 1:length(good_colour)
            this_colour = good_colour(colour_ix);
            for ix = 1:2
                this_bool = corner_bool(:, ix, this_colour);
                firingRates(:, ix, colour_ix, :, 1:sum(this_bool)) = newD(:, :, this_bool);
            end
        end
        
        if length(good_colour) > 1
            combinedParams = {{[1 3]}, {1}, {2}, {3}};
            margNames = {'Context-Time', 'Context', 'Colour', 'Time'};
        else
            firingRates = squeeze(firingRates);
            combinedParams = {{[1 2]}, {1}, {2}};
            margNames = {'Context-Time', 'Context', 'Time'};
        end
        
        firingRatesAverage = nanmean(firingRates, ndims(firingRates));
        
        %optimalLambda = 0.0025;
        [W,V,whichMarg] = dpca(firingRatesAverage, min([nUnits nTimes 15]), ...
                                'combinedParams', combinedParams);
                                %'lambda', optimalLambda);
        
        %Uncomment the following to plot dPCA results.
        explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
            'combinedParams', combinedParams);
        dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
            'explainedVar', explVar, ...
            'marginalizationNames', margNames, ...
            'whichMarg', whichMarg,                 ...
            'timeMarginalization', 2 + double(length(good_colour) > 1), ...
            'legendSubplot', 12);
        
    end
    
    
    
end
