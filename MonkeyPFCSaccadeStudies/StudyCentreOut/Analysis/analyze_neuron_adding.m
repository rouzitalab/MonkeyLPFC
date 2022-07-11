addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths DEF
addpath(genpath(fullfile(paths.ml3rd, 'DataHigh')));
addpath(fullfile(paths.ml3rd, 'dPCA', 'matlab'));

%%
% %% BCILAB - startup is very slow.
% %TODO: Get rid of this once I find better rLDA and rLogReg.
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

%% Constants
sessions = my_sessions('CentreOut');
analysisParams = my_anaparams('classification');
anaWin = analysisParams.anaWins(strcmpi({analysisParams.anaWins.name}, 'fullTrial'));
analysisParams.anaWins = anaWin;
nAna = length(analysisParams.anaWins);

tuningParams = my_anaparams('tuning');  % Will be used for zerovar epochs

doStandardizeX = true;

%% Prepare output variables
neur_order = cell(1, length(sessions));
neur_acc = cell(1, length(sessions));
rneur_acc = cell(1, length(sessions));
nReps = 20;

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
    
    %% Prepare cross-validations so classes are well distributed across folds.
    % We will use these same CVs for feature reduction and for machine
    % learning.
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
    
    %% Prepare data
    trainX = cell(analysisParams.kcross, 1);
    testX = trainX;
    
    cDat = trimTrials(cDat, anaWin);
    D = cDat2DataHigh(cDat.trial);
    
    %Remove neurons with one or more bad epochs
    for tr_ix = 1:length(D)
        D(tr_ix).data = D(tr_ix).data(keepBool, :);
    end
    
    % Convert to firing rate from entire trial.
    for tr_ix = 1:length(D)
        D(tr_ix).data = 1000 * nansum(D(tr_ix).data, 2) ./ sum(~isnan(D(tr_ix).data), 2);
    end
    
    % Reset the DataHigh parameters
    dhParams = struct(...
        'binWidth', 1, ...%analysisParams.binWidth,...
        'use_sqrt', false, ...%true,...
        'kern', 0, ...%analysisParams.kernSD,...
        'trial_average', false,...
        'keep_neurons', cDat.invalidUnits(keepBool) == DEF.unitState.valid);
    nDimRed = sum(dhParams.keep_neurons);
    dhRedMode = -1;
    
    %Get the (reduced) feature set
    [newD, ~, ~] = reducedims(D, dhRedMode, nDimRed, dhParams);
    
    %Cut to equal length trials.
    dataL = min(cellfun('size', {newD.data}, 2));
    for tr_ix = 1:length(newD)
        newD(tr_ix).data = newD(tr_ix).data(:, 1:dataL);
    end; clear tr_ix
    
    %-> trials x neurons x samples
    newD = permute(cat(3, newD.data), [3 1 2]);
    [nTrials, nNeurons, nSamples] = size(newD);
    
    % Identify neurons that may have zero variance in one of the CVs
    keepFeat = true(1, nNeurons);
    for cv_ix = 1:analysisParams.kcross
        trainBool = ~testBool(:, cv_ix);
        
        thisTrainX = newD(trainBool, :);
        thisTestX = newD(~trainBool, :);
        
        keepFeat = keepFeat & std([thisTrainX;thisTestX]) > 0 & ...
            var(thisTrainX) > 0;
    end
    
    newD = newD(:, keepFeat);
    
    if doStandardizeX
        newD = zscore(newD);
    end
    
    % Put data into cell arrays for train and test.
    for cv_ix = 1:analysisParams.kcross
        trainBool = ~testBool(:, cv_ix);
        
        trainX{cv_ix} = newD(trainBool, :);
        testX{cv_ix} = newD(~trainBool, :);
        
    end; clear cv_ix
    
%     %% Greedy neuron-adding classification accuracy
% %     addpath(bcipath{:});
%     
%     neur_used = [];
%     neur_avail = 1:nNeurons;
%     acc_out = nan(1, nNeurons);
%     
%     for loop_ix = 1:nNeurons
%         
%         % Test the addition of each remaining neuron to see which contributes
%         % most to cross-validated classification accuracy.
%         loop_acc = nan(1, length(neur_avail));
%         for neur_ix = 1:length(neur_avail)
%             neur_ids = [neur_used neur_avail(neur_ix)];
%             nNeurs = length(neur_ids);
%             %cross-validated
%             this_predY = nan(size(Ynum));
%             for cv_ix = 1:analysisParams.kcross
%                 trainBool = ~testBool(:, cv_ix);
%                 
%                 nTrain = sum(trainBool);
%                 nTest = sum(~trainBool);
%                 
%                 thisTrainY = Ynum(trainBool);
%                 % Get training data and reshape so each time-point is a
%                 % different variable.
%                 thisTrainX = reshape(trainX{cv_ix}(:, neur_ids, :), [nTrain, nNeurs*nSamples]);
%                 thisTestX = reshape(testX{cv_ix}(:, neur_ids, :), [nTest, nNeurs*nSamples]);
%                 
%                 % Do the machine learning.
%                 if ~isempty(thisTrainX) && ~isempty(thisTestX) && ~isempty(thisTrainY)
%                     
%                     %BCILAB method.
% %                     model = ml_train({[ones(size(thisTrainX, 1), 1) thisTrainX], thisTrainY},...
% %                         {'lda',...
% %                         'lambda', [],...
% %                         'regularization', 'auto'});  %'auto' 'shrinkage' 'independence'
% %                     predictions = ml_predict([ones(size(thisTestX, 1), 1) thisTestX], model);
% %                     [~, this_predY(~trainBool)] = max(predictions{2}, [], 2);
%                     
%                     linclass = fitcdiscr(thisTrainX, thisTrainY);
%                     this_predY(~trainBool) = predict(linclass, thisTestX);
%                 end
%                 clear nTrain nTest thisTrainY thisTrainX thisTestX model trainBool
%             end; clear cv_ix
%             loop_acc(neur_ix) = 100*sum(this_predY == Ynum)./nTrials;
%         end; clear neur_ix
%         [acc_out(loop_ix), best_neur_ix] = max(loop_acc);
%         neur_used = cat(2, neur_used, neur_avail(best_neur_ix));
%         neur_avail(best_neur_ix) = [];
%     end; clear loop_ix
%     
%     neur_order{sess_ix} = neur_used;
%     neur_acc{sess_ix} = acc_out;
%     
% %     rmpath(bcipath{:});

    %% Random neuron adding
    acc_out = nan(nReps, nNeurons);
    for rep_ix = 1:nReps
        neur_used = false(1, nNeurons);
        for neur_ix = 1:nNeurons
            neur_avail = find(~neur_used);
            this_neur_id = neur_avail(randperm(length(neur_avail), 1));
            neur_used(this_neur_id) = true;
            
            %cross-validated
            this_predY = nan(size(Ynum));
            for cv_ix = 1:analysisParams.kcross
                trainBool = ~testBool(:, cv_ix);
                
                nTrain = sum(trainBool);
                nTest = sum(~trainBool);
                
                thisTrainY = Ynum(trainBool);
                % Get training data and reshape so each time-point is a
                % different variable.
                thisTrainX = reshape(trainX{cv_ix}(:, neur_used, :), [nTrain, sum(neur_used)*nSamples]);
                thisTestX = reshape(testX{cv_ix}(:, neur_used, :), [nTest, sum(neur_used)*nSamples]);
                
                % Do the machine learning.
                if ~isempty(thisTrainX) && ~isempty(thisTestX) && ~isempty(thisTrainY)
                    
                    %BCILAB method.
%                     model = ml_train({[ones(size(thisTrainX, 1), 1) thisTrainX], thisTrainY},...
%                         {'lda',...
%                         'lambda', [],...
%                         'regularization', 'auto'});  %'auto' 'shrinkage' 'independence'
%                     predictions = ml_predict([ones(size(thisTestX, 1), 1) thisTestX], model);
%                     [~, this_predY(~trainBool)] = max(predictions{2}, [], 2);
                    
                    linclass = fitcdiscr(thisTrainX, thisTrainY);
                    this_predY(~trainBool) = predict(linclass, thisTestX);
                end
                clear nTrain nTest thisTrainY thisTrainX thisTestX model trainBool
            end; clear cv_ix
            acc_out(rep_ix, neur_ix) = 100*sum(this_predY == Ynum)./nTrials;
        end; clear neur_ix 
    end; clear rep_ix
    rneur_acc{sess_ix} = acc_out;
    
    %%
    fprintf(['Finished ' num2str(sess_ix) ' of ' num2str(length(sessions)) '.\n']);
    clear testBool ptb cDat trialBool Ynum Ystr nTrials this_sess
    clear neur_used neur_avail acc_out
end
save(fullfile(paths.results, 'neuron_adding.mat'), ...
    'neur_order', 'neur_acc', 'rneur_acc', 'analysisParams', 'sessions');