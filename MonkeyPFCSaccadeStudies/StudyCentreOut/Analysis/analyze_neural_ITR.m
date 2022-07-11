%% Analyze information transfer rate
% using only neural trajectories and rLDA
% at different trial lengths.

addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths DEF
addpath(genpath(fullfile(paths.ml3rd, 'DataHigh')));
addpath(fullfile(paths.ml3rd, 'dPCA', 'matlab'));

%% BCILAB - startup is very slow.
% prevwd = pwd();
% bcipath = fullfile(paths.ml3rd,'BCILAB');
% addpath(bcipath);
% bcilab;
% cd(prevwd);
% p = strsplit(path, pathsep);
% isbcipath = false(size(p));
% for p_ix = 1:length(p)
%     isbcipath(p_ix) = length(p{p_ix}) >= length(bcipath) && strcmpi(p{p_ix}(1:length(bcipath)), bcipath);
% end
% bcipath = p(isbcipath);
% rmpath(bcipath{:});
% clear prevwd p isbcipath p_ix CURRENTSTUDY

%% Paths & Constants
sessions = my_sessions('CentreOut');
analysisParams = my_anaparams('ITR'); % Only trajectory
nAna = length(analysisParams.anaWins);  % Should be == 1
anaWin = analysisParams.anaWins;

tuningParams = my_anaparams('tuning');  % Will be used for zerovar epochs

% Feature reduction parameters (using datahigh)
dhParams = struct(...
    'binWidth', analysisParams.binWidth,...
    'use_sqrt', true,...
    'kern', analysisParams.kernSD,...
    'trial_average', false,...
    'keep_neurons', true(1,1));

% Feature set names. Use 'smooth' instead of 'trajectory'.
featureNames = {'fullTrial'};
nFeatureSets = length(featureNames);

stepSize = dhParams.binWidth;  % msec
winEdges = analysisParams.anaWins(1).winEdges;
winStops = winEdges(1)+stepSize:stepSize:winEdges(end)+stepSize;
nWinSteps = length(winStops);

doStandardizeX = true;

%% Prepare output variables
predictedY = cell(length(sessions), 1);
actualY = cell(length(sessions), 1);
actualYstr = actualY;

sess_ix = 1;  % Quick init for debug.
%%
for sess_ix = 1:length(sessions)
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
    
    % Eliminate units that fire < 1Hz
    % Adds cDat.invalidUnits.
    % If thisParams.useUnits == 'merged' then it merges all units per
    % channel into one unit per channel.
    cDat = triageUnits(cDat, analysisParams);
    
    %% Dependent Variables (Y)
    % [nTrials x 1]
    Ynum = cat(1, cDat.trial.newClass);
    Ystr = {cDat.trial.newClassStr}';
    nTrials = length(Ynum);
    
    actualY{sess_ix} = Ynum;
    actualYstr{sess_ix} = Ystr;
    
    %% Prepare cross-validations so classes are well distributed across folds.
    testBool = cvDistrib(Ynum, analysisParams.kcross);
    
    %% Identify neurons with 0 variance across trials in any epoch as 'bad'
    % Do for each CV
    nEpochs = length(tuningParams.anaWins);
    nNeurons = length(cDat.invalidUnits);
    goodEpochCVBool = false(nNeurons, nEpochs, size(testBool, 2));
    for ep_ix = 1:nEpochs
        ana_win = tuningParams.anaWins(ep_ix);
        raster = getEpochRaster(cDat, ana_win);  % trials x units x samples
        nGoodSamps = sum(~isnan(raster), 3);
        nSpikes = sum(raster > 0, 3);
        frate = 1000*nSpikes./nGoodSamps;
        for cv_ix = 1:size(testBool, 2)
            trainBool = ~testBool(:, cv_ix);
            goodEpochCVBool(:, ep_ix, cv_ix) = var(frate(trainBool, :)) > 0;
        end
    end
    keepBool = ~any(any(~goodEpochCVBool, 3), 2);
    clear goodEpochCVBool ep_ix ana_win raster nGoodSamps nSpikes frate cv_ix
    
    %% For each window size
    for win_ix = 1:nWinSteps
        
        %Trim trials for this window size.
        tmp_ana_win = anaWin;
        tmp_ana_win.winEdges = [winEdges(1) winStops(win_ix)];
        cDat = trimTrials(cDat, tmp_ana_win);
        clear tmp_ana_win
        
        % Get the data in DataHigh format
        D = cDat2DataHigh(cDat.trial);
        
        %Remove neurons with one or more bad epochs
        for tr_ix = 1:length(D)
            D(tr_ix).data = D(tr_ix).data(keepBool, :);
        end
        
        % Reset the DataHigh parameters
        dhParams.keep_neurons = cDat.invalidUnits(keepBool) == DEF.unitState.valid;
        nDimRed = sum(dhParams.keep_neurons);
        dhRedMode = -1; % Just smoothing
        
        tempParams = dhParams;
        if strcmpi(featureNames, 'fullTrial')
            for tr_ix = 1:length(D)
                D(tr_ix).data = 1000 * nansum(D(tr_ix).data, 2) ./ sum(~isnan(D(tr_ix).data), 2);
            end
            tempParams.use_sqrt = false;
            tempParams.binWidth = 1;
            tempParams.kern = 0;
        else
            
        end
        [redD, C, lat, dhRedParams] = reducedims(D, dhRedMode, nDimRed, tempParams);
        
        % Trim all trials
        dataL = min(cellfun('size', {redD.data}, 2));
        for tr_ix = 1:length(redD)
            redD(tr_ix).data = redD(tr_ix).data(:, 1:dataL);
        end; clear tr_ix
        
        predY = nan(size(Ynum));
%         addpath(bcipath{:});
        for cv_ix = 1:analysisParams.kcross
            trainBool = ~testBool(:, cv_ix);
            nTrain = sum(trainBool);
            nTest = sum(~trainBool);
            
            thisTrainY = Ynum(trainBool);
            thisTrainX = permute(cat(3, redD(trainBool).data), [3 1 2]);
            thisTestX = permute(cat(3, redD(~trainBool).data), [3 1 2]);
            [~, nNeurons, nTimes] = size(thisTrainX);
            thisTrainX = reshape(thisTrainX, [nTrain, nNeurons*nTimes]);
            thisTestX = reshape(thisTestX, [nTest, nNeurons*nTimes]);
            
            if doStandardizeX
                x_bar = mean([thisTrainX;thisTestX]);
                x_std = std([thisTrainX;thisTestX]);

                keepFeat = x_std > 0 & var(thisTrainX) > 0;

                thisTrainX = bsxfun(@rdivide, bsxfun(@minus, thisTrainX, x_bar),x_std);
                thisTestX = bsxfun(@rdivide, bsxfun(@minus, thisTestX, x_bar), x_std);


                thisTrainX = thisTrainX(:, keepFeat);
                thisTestX = thisTestX(:, keepFeat);
            end
            
            %% Classification
%             model = ml_train({[ones(size(thisTrainX, 1), 1) thisTrainX], thisTrainY},...
%                 {'lda',...
%                 'lambda', [],...
%                 'regularization', 'auto'});  %'auto' 'shrinkage' 'independence'
%             predictions = ml_predict([ones(size(thisTestX, 1), 1) thisTestX], model);
%             [~, predY(~trainBool)] = max(predictions{2}, [], 2);
            
            linclass = fitcdiscr(thisTrainX, thisTrainY);%, 'discrimType', 'quadratic');
            predY(~trainBool) = predict(linclass, thisTestX);
            
        end; clear cv_ix  % for cv_ix
%         rmpath(bcipath{:});
        
        predictedY{sess_ix}(:, win_ix) = predY;
        clear predY
        
    end; clear win_ix

    
    fprintf(['Finished ' num2str(sess_ix) ' of ' num2str(length(sessions)) '.\n']);
end
save(fullfile(paths.results, 'ITR.mat'), ...
    'actualY', 'actualYstr', 'analysisParams', 'dhParams',...
    'predictedY', 'sessions', 'winStops');