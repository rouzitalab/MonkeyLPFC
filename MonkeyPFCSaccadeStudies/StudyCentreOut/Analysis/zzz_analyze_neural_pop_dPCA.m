%% Paths and Constants
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths
sessions = my_sessions('CentreOut');
analysisParams = my_anaparams('dPCA');

%addpath(fullfile(paths.mlDevel, 'Helper'));
addpath(genpath(fullfile(paths.ml3rd, 'DataHigh')));
addpath(fullfile(paths.ml3rd, 'dPCA', 'matlab'));

%% Constants and other frequently accessed variables
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
trial_frate = cell(nSessions, nAna);  % each cell is nNeurons x nTimes x nClasses x nTrials
tVec = cell(nSessions, nAna);

%nTimes and nClasses will have to be permuted prior to dPCA

sess_ix = 2;  % Need to replace session 2 because its class==1 trials have almost no post-saccade data.
%%
for sess_ix = 1:nSessions
    %%
    this_sess = sessions(sess_ix);
    
    %% Load and process ptbmat
    ptb = load(fullfile(paths.preprocessed, 'ptb', this_sess.ptbfname), '-mat');
    [ptb.eData.trial, ptb.params] = getTrialStimInfo(ptb.eData.trial, ptb.params);
    ptb.eData.trial = getTrialBehavResult(ptb.eData.trial, ptb.params, analysisParams);
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
    
    % Eliminate units that fire < 1Hz or > 100 Hz.
    % Adds cDat.invalidUnits.
    % If thisParams.useUnits == 'merged' then it merges all units per
    % channel into one unit per channel.
    cDat = triageUnits(cDat, analysisParams);
    
    %% Easy-access, map classes -> directions
    classId = [cDat.trial.newClass];
    [uqClasses, inClassIx] = unique(classId);
    classCount = hist(classId, uqClasses);
    nClasses = length(uqClasses);
    clear inClassIx
    
    trial_class{sess_ix, 1} = classId';
    
    nTrials = length(cDat.trial);
    nUnits = size(cDat.trial(1).raster, 2);
    dhParams.keep_neurons = true(nUnits, 1);
    sess_unit_id{sess_ix} = [sess_ix*ones(nUnits,1) (1:nUnits)'];
    
    trialNum{sess_ix, 1} = ones(nUnits, 1) * classCount;
    
    %%
    for win_ix = 1:nAna
        
        %%
        ana_win = analysisParams.anaWins(win_ix);
        
        % Identify eventTime, avoidTime, a tVec, and timeBool for each trial.
        cDat = trimTrials(cDat, ana_win);
        
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
        
        %Save the tVec as an output
        tr_ix = find(dSize == max(dSize), 1, 'first');
        tempT = cDat.trial(tr_ix).tVec(cDat.trial(tr_ix).timeBool);
        tempT = tempT(1, 1:dhParams.binWidth * floor(size(tempT, 2)/dhParams.binWidth));  % Cut off excess
        tempT = reshape(tempT, dhParams.binWidth, []);
        tempT = mean(tempT);
        tVec{sess_ix, win_ix} = tempT;
        
        %datahigh2dpca(newD)
        newD = cat(3, newD.data);
        nTimes = size(newD, 2);
        
        % Allocate memory
        trial_frate{sess_ix, win_ix} = nan(nUnits, nTimes, nClasses, max(classCount));
        
        for cl_ix = 1:nClasses
            trBool = classId == cl_ix;
            trial_frate{sess_ix, win_ix}(:, :, cl_ix, 1:sum(trBool)) = newD(:, :, trBool);
        end
    end
    
end
trialNum = cat(1, trialNum{:});
sess_unit_id = cat(1, sess_unit_id{:});
trial_class = cat(1, trial_class{:});

%% Prepare for dPCA
mySizes = cell2mat(cellfun(@size, trial_frate, 'UniformOutput', false));
maxTrialNum = max(mySizes(:, end));  % Max number of trials.
N = sum(mySizes(:,1));  % Number of Neurons
S = nClasses;  % Number of Stimuli
T = sum(mySizes(1, [2 6]));  % Number of Time steps.
firingRates = nan(N, S, T, maxTrialNum);
for sess_ix = 1:nSessions
    temp = cat(2, trial_frate{sess_ix, :});
    unitBool = sess_unit_id(:, 1) == sess_ix;
    firingRates(unitBool, :, :, 1:size(temp, 4)) = permute(temp, [1 3 2 4]);
end
firingRatesAverage = nanmean(firingRates, 4);  % Average across trials.
nTrials = sum(~isnan(firingRates), 4);
clear mySizes temp sess_ix unitBool

% Cut out any time points that are missing data for all trials for any
% given neuron-location combination.
goodTimesBool = squeeze(sum(sum(isnan(firingRatesAverage), 2), 1) == 0);
firingRatesAverage = firingRatesAverage(:, :, goodTimesBool);
firingRates = firingRates(:, :, goodTimesBool, :);

%% Setup dPCA

combinedParams = {{2}, {1, [1 2]}};
margNames = {'Time', 'Location'};

% The following generates the time-only, and time*location interaction
% components. The latter show trials starting at different offsets then
% converging toward the middle then diverging again.
% combinedParams = {{2}, {[1 2]}};
% margNames = {'Time', 'Location*Time'};

% The following generates time-only, location-only (horizontal
% lines at different offsets), and time-location interaction components.
% The latter are annoying because they mostly undo the offsets generated by
% the location-only components, only at different times.
% combinedParams = {{2}, {1}, {[1 2]}};
% margNames = {'Time', 'Location', 'Location*Time'};


%% dPCA
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization

% optimalLambda = 0.0025;
%uncomment the following to calculate optimalLambda
%(a bit slow).
% optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
%     'combinedParams', combinedParams, ...
%     'numComps', 20,...
%     'numRep', 10, ...  % increase this number to ~10 for better accuracy
%     'filename', 'tmp_optimalLambdas.mat');

[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams, 'scale', 'yes');%,...
%     'lambda',optimalLambda);

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

beep;
save(fullfile(paths.results, 'dPCA_dat.mat'), ...
    'analysisParams', 'nTrials', 'tVec', 'sess_unit_id',...
    'firingRatesAverage', 'combinedParams', 'margNames', 'W', 'V', 'whichMarg', 'explVar');

%% From tutorial
% %% Step 1: PCA of the dataset
% X = firingRatesAverage(:,:);
% X = bsxfun(@minus, X, mean(X,2));
% [W,~,~] = svd(X);
% 
% time = tVec{1, 1};
% time = [time time(end)+50+abs(tVec{1, 2}(1))+tVec{1,2}];
% timeEvents = [0 0; 250 250; 1300 1300];
% 
% % minimal plotting
% % dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...);
% 
% % computing explained variance
% explVar = dpca_explainedVariance(firingRatesAverage, W, W, ...
%     'combinedParams', combinedParams);
% 
% % a bit more informative plotting
% margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
% dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...
%     'explainedVar', explVar, ...
%     'marginalizationNames', margNames, ...
%     'marginalizationColours', margColours, ...
%     'time', time, ...
%     'timeEvents', timeEvents);
% 
% 
% %% Step 2: PCA in each marginalization separately
% dpca_perMarginalization(firingRatesAverage, @dpca_plot_default, ...
%    'combinedParams', combinedParams);
% 
% 
% %% Step 3: dPCA without regularization
% 
% % This is the core function.
% % W is the decoder, V is the encoder (ordered by explained variance),
% % whichMarg is an array that tells you which component comes from which
% % marginalization
% 
% [W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
%     'combinedParams', combinedParams);
% 
% explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
%     'combinedParams', combinedParams);
% 
% dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
%     'explainedVar', explVar, ...
%     'marginalizationNames', margNames, ...
%     'marginalizationColours', margColours, ...
%     'whichMarg', whichMarg,                 ...
%     'time', time,                        ...
%     'timeEvents', timeEvents,               ...
%     'timeMarginalization', 2, ...
%     'legendSubplot', 12);


% %% Step 4: dPCA with regularization
% 
% % This function takes some minutes to run. It will save the computations 
% % in a .mat file with a given name. Once computed, you can simply load 
% % lambdas out of this file:
% %   load('tmp_optimalLambdas.mat', 'optimalLambda')
% 
% optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
%     'combinedParams', combinedParams, ...
%     'numRep', 5, ...  % increase this number to ~10 for better accuracy
%     'filename', 'tmp_optimalLambdas.mat');
% 
% [W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
%     'combinedParams', combinedParams, ...
%     'lambda', optimalLambda);
% 
% explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
%     'combinedParams', combinedParams);
% 
% dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
%     'explainedVar', explVar, ...
%     'marginalizationNames', margNames, ...
%     'marginalizationColours', margColours, ...
%     'whichMarg', whichMarg,                 ...
%     'time', time,                        ...
%     'timeEvents', timeEvents,               ...
%     'timeMarginalization', 3,           ...
%     'legendSubplot', 16);
% 
% 
% % sess_unit_id([84 103 78 58 85], :)
% 
% %% Decoding
% 
% decodingClasses = {[(1:S)' (1:S)'], repmat([1:2], [S 1]), [], [(1:S)' (S+(1:S))']};
% 
% accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
%     'lambda', optimalLambda, ...
%     'combinedParams', combinedParams, ...
%     'decodingClasses', decodingClasses, ...
%     'numRep', 5, ...        % increase to 100
%     'filename', 'tmp_classification_accuracy.mat');
% 
% accuracyShuffle = dpca_classificationShuffled(firingRates, trialNum, ...
%     'lambda', optimalLambda, ...
%     'combinedParams', combinedParams, ...
%     'decodingClasses', decodingClasses, ...
%     'numRep', 5, ...        % increase to 100
%     'numShuffles', 20, ...  % increase to 100 (takes a lot of time)
%     'filename', 'tmp_classification_accuracy.mat');
% 
% componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);
% 
% dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
%     'explainedVar', explVar, ...
%     'marginalizationNames', margNames, ...
%     'marginalizationColours', margColours, ...
%     'whichMarg', whichMarg,                 ...
%     'time', time,                        ...
%     'timeEvents', timeEvents,               ...
%     'timeMarginalization', 3,           ...
%     'legendSubplot', 16,                ...
%     'componentsSignif', componentsSignif);