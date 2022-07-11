%%
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths DEF
addpath(genpath(fullfile(paths.ml3rd, 'DataHigh')));
% addpath(fullfile(paths.ml3rd, 'dPCA', 'matlab'));

%% BCILAB - startup is very slow.
% % %TODO: Get rid of this once I find better rLDA and rLogReg.
% prevwd = pwd();
% bcipath = fullfile(paths.ml3rd, 'BCILAB');
% addpath(bcipath);
% bcilab;
% cd(prevwd);
% p = strsplit(path, ':');
% isbcipath = false(size(p));
% for p_ix = 1:length(p)
%     isbcipath(p_ix) = length(p{p_ix}) >= length(bcipath) && strcmpi(p{p_ix}(1:length(bcipath)), bcipath);
% end
% bcipath = p(isbcipath);
% rmpath(bcipath{:});
% clear prevwd p isbcipath p_ix

%% Constants
sessions = my_sessions('RegionRule');
analysisParams = my_anaparams('classification');
nAna = length(analysisParams.anaWins);  % Number of analysis windows.
% analysisParams.anaWins(strcmpi({analysisParams.anaWins.name}, 'trajectory')) = [];

tuningParams = my_anaparams('tuning');  % Will be used for zerovar epochs

% Feature reduction parameters (using datahigh)
dhParams = struct(...
    'binWidth', analysisParams.binWidth,...
    'use_sqrt', true,...
    'kern', analysisParams.kernSD,...
    'trial_average', false,...
    'keep_neurons', true(1,1));

% Feature set names
featureNames = {analysisParams.anaWins.name};
%'baseline'    'target'    'cue'    'delay'    'fullDelay'    'response'    'fullTrial'    'trajectory'
trajNames = {'smooth' 'gpfa'};
traj_bool = strcmpi(featureNames, 'trajectory');
if any(traj_bool)
    featureNames(traj_bool) = [];
    featureNames = [featureNames trajNames];
end
clear traj_bool

featureNames = [featureNames 'baseline,cue,delay,response']; %Concatenated spike counts.

nFeatureSets = length(featureNames);

% Machine learning names.
mlNames = {'SVM-Lin'};  %'SVM-Lin', 'SVM-RBF', 'LDA', 'rLDA'};  %DBN, DTW, NB
nML = length(mlNames);

maxTrajDims = 7;
nTargPairs = 4;
doStandardizeX = true;
d_thresh = 1.5;

%% Prepare output variables
predictedY = cell(length(sessions), nFeatureSets, nML);
actualY = cell(length(sessions), 1);
actualYstr = actualY;
pairBoolOut = cell(length(sessions), 1);

sess_ix = 1;
%%
for sess_ix = 1:length(sessions)
    %%
    this_sess = sessions(sess_ix);
    
    %% Load and process ptbmat
    ptb = load(fullfile(paths.preprocessed, 'ptb', this_sess.ptbfname), '-mat');
    [ptb.eData.trial, ptb.params] = getTrialStimInfo(ptb.eData.trial, ptb.params);
    ptb.eData.trial = getTrialBehavResult(ptb.eData.trial, ptb.params, analysisParams);
    [IsInCondition] = getCornerConditionMembership(ptb.eData.trial);
    ptb.eData.trial = getNewClass(ptb.eData.trial, analysisParams);
    behav = procBehavBlocks(ptb.eData.trial);
    pmarg = 0.01;
    tpr = behav.TrueA./(behav.TrueA + behav.FalseB);
    tpr(tpr<pmarg)=pmarg; tpr(tpr>1-pmarg)=1-pmarg;
    fpr = behav.FalseA./(behav.FalseA+behav.TrueB);
    fpr(fpr<pmarg)=pmarg; fpr(fpr>1-pmarg)=1-pmarg;
    ztpr = norminv(tpr, 0, 1);
    zfpr = norminv(fpr, 0, 1);
    d = ztpr - zfpr;
    clear behav pmarg tpr fpr ztpr zfpr
    
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
    IsInCondition = IsInCondition(trialBool, :, :, :);
    d = d(trialBool);
    
    % Eliminate units that...
    % Adds cDat.invalidUnits.
    % If thisParams.useUnits == 'merged' then it merges all units per
    % channel into one unit per channel.
    cDat = triageUnits(cDat, analysisParams);
    
    
     %% Dependent Variables (Y)
    % [nTrials x 1]
    trClass = cat(1, cDat.trial.newClass);
    trStr = {cDat.trial.newClassStr}';
    
    %% Find matching corners
    pairBool = false(sum(trialBool), 2, 4);
    for corner_ix = 1:4
        pairBool(:, :, corner_ix) = squeeze(any(any(IsInCondition(:, corner_ix, :, :, :), 5), 4));
    end
    has_pairing = ~any(~(squeeze(sum(pairBool)) > 30));
    pairBool = pairBool(:, :, has_pairing);
    pairBoolOut{sess_ix} = pairBool;
    
    %% For each pairing
    for pair_ix = 1:size(pairBool, 3)
        thisBool = pairBool(:, :, pair_ix);
        thisBool = bsxfun(@and, thisBool, d > d_thresh);
        
        eitherBool = any(thisBool, 2);
        smallTestBool = cvDistrib(trClass(eitherBool), analysisParams.kcross);
        testBool = false(length(cDat.trial), analysisParams.kcross);
        trainBool = testBool;
        testBool(eitherBool, :) = smallTestBool;
        trainBool(eitherBool, :) = ~smallTestBool;
        
        actualY{sess_ix}{pair_ix} = trClass(eitherBool);
        actualYstr{sess_ix}{pair_ix} = trStr(eitherBool);
        
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
                goodEpochCVBool(:, ep_ix, cv_ix) = var(frate(trainBool(:, cv_ix), :)) > 0;
            end
        end
        goodNeuronBool = ~any(any(~goodEpochCVBool, 3), 2);
        clear goodEpochCVBool ep_ix ana_win raster nGoodSamps nSpikes frate cv_ix
        
        %% For each feature set
        for fs_ix = 1:nFeatureSets
            
            fsName = featureNames{fs_ix};
            isTraj = any(strcmpi(trajNames, fsName));
            subfsNames = strsplit(fsName, ',');
            
            subTrainX = cell(1, length(subfsNames));
            subTestX = cell(1, length(subfsNames));
            
            for sfs_ix = 1:length(subfsNames)
                sfsNames = strsplit(subfsNames{sfs_ix}, '_'); %split again in case canonical
                sfsName = sfsNames{1};
                
                % Convert the data for dataHigh
                if isTraj
                    ana_win = analysisParams.anaWins(strcmpi({analysisParams.anaWins.name}, 'trajectory'));
                else
                    ana_win = analysisParams.anaWins(strcmpi({analysisParams.anaWins.name}, sfsName));
                end
                cDat = trimTrials(cDat, ana_win);
                timeVec = cDat.trial(1).tVec(cDat.trial(1).timeBool);
                
                D = cDat2DataHigh(cDat.trial);
                clear ana_win
                
                %Remove neurons with one or more bad epochs
                for tr_ix = 1:length(D)
                    D(tr_ix).data = D(tr_ix).data(goodNeuronBool, :);
                end
                
                % Reset the DataHigh parameters
                dhParams.binWidth = analysisParams.binWidth;
                dhParams.kern = analysisParams.kernSD;
                dhParams.keep_neurons = cDat.invalidUnits(goodNeuronBool) == DEF.unitState.valid;
                nDimRed = sum(dhParams.keep_neurons);
                dhRedMode = -1;
                
                trainX = cell(analysisParams.kcross, 1);
                testX = trainX;
                
                %TODO: Standardize data by baseline PER BLOCK
                if ~isTraj
                    
                    if strcmpi(sfsName, 'fullTrial')
                        %Get the average firing rate for each trial
                        for tr_ix = 1:length(D)
                            D(tr_ix).data = 1000 * nansum(D(tr_ix).data, 2) ./ sum(~isnan(D(tr_ix).data), 2);
                        end
                        tempParams = dhParams;
                        tempParams.use_sqrt = false;
                        tempParams.binWidth = 1;
                        tempParams.kern = 0;
                    else
                        % Get the sqrt spike count for the epoch
                        lcdLength = nan(size(D));
                        for tr_ix = 1:length(D)
                            lcdLength(tr_ix) = size(D(tr_ix).data, 2);
                        end; clear tr_ix
                        tempParams = dhParams;
                        tempParams.binWidth = min(lcdLength);
                        tempParams.kern = 0;
                        clear lcsdLength
                    end
                    [newD, ~, ~] = reducedims(D, dhRedMode, nDimRed, tempParams);
                    %Cut to equal length trials.
                    dataL = min(cellfun('size', {newD.data}, 2));
                    for tr_ix = 1:length(newD)
                        newD(tr_ix).data = newD(tr_ix).data(:, 1:dataL);
                    end; clear tr_ix
                    
                    %-> trials x neurons x samples
                    newD = permute(cat(3, newD.data), [3 1 2]);
                    [nTrials, nNeurons, nSamples] = size(newD);
                    for cv_ix = 1:analysisParams.kcross
                        trainX{cv_ix} = newD(trainBool(:, cv_ix), :, :);
                        testX{cv_ix} = newD(testBool(:, cv_ix), :, :);
                    end; clear cv_ix
                
                else
                    %set dhRedMode
                    %smooth: -1; pca:1; fa:3; lda: 4, gpfa:5
                    %Default is smooth
                    if strcmpi(sfsName, 'pca')
                        dhRedMode = 1;
                        nDimRed = maxTrajDims;
                    elseif strcmpi(fsName, 'fa')
                        dhRedMode = 3;
                        nDimRed = maxTrajDims;
                    elseif strcmpi(fsName, 'lda')
                        dhRedMode = 4;
                        nDimRed = length(unique({D.condition}));
                    elseif strcmpi(sfsName, 'gpfa')
                        dhRedMode = 5;
                        dhParams.kern = 0;
                        nDimRed = maxTrajDims;
                    end
                    
                    % TODO: Standardize by baseline per block
                    % but (GP)FA does not accept mean-0 data.
                    

                    % Trajectory transformations must be calculated using training
                    % data only.
                    for cv_ix = 1:analysisParams.kcross
                        
                        [trainD, C, lat, dhRedParams] =...
                            reducedims(D(trainBool(:, cv_ix)), dhRedMode, nDimRed, dhParams);
                        
                        % Then test data are transformed using learned
                        % parameters
                        %Convert testing dat to transformable format
                        [testD, ~, ~, ~] =...
                            reducedims(D(testBool(:, cv_ix)), -1, size(D(1).data, 1), dhParams);
                        [testD.y] = deal(testD.data);  % Copy data into y.
                        t_length = nan(size(testD));
                        for tr_ix = 1:numel(testD)
                            t_length(tr_ix) = size(testD(tr_ix).data, 2);
                        end
                        testIds = find(testBool(:, cv_ix));
                        tId = num2cell(testIds);
                        [testD.trialId] = tId{:};
                        t_length = num2cell(t_length);
                        [testD.T] = t_length{:};
                        clear tr_ix t_length
                    
                        % Transform testD.data
                        if any(strcmpi({'pca', 'fa', 'lda'}, sfsName))
                            offsets = mean([testD.data], 2);
                            %offsets = zeros(size([testD.data], 1), 1);
                            for tr_ix = 1:length(testD)
                                testD(tr_ix).data = C' * bsxfun(@minus, testD(tr_ix).data, offsets);
                            end; clear tr_ix offsets
                        elseif strcmpi(sfsName, 'gpfa')
                            testD = exactInferenceWithLL(testD, dhRedParams);
                            [Xorth, Corth] = orthogonalize([testD.xsm], dhRedParams.C);
                            testD = segmentByTrial(testD, Xorth, 'data');
                            testD = rmfield(testD, {'Vsm', 'VsmGP', 'xsm'});
                            clear Xort Corth
                        end
                        clear testIds tId s C lat
                        
                        %% Convert structure to output matrices
                        
                        % Reformat epochStarts to match the binWidth.
                        % Unnecessary because it is unused.
                        %                 for tr_ix = 1:length(newD)
                        %                     newD(tr_ix).epochStarts = ceil(newD(tr_ix).epochStarts / dhParams.binWidth);
                        %                 end
                        
                        % Trim all trials
                        dataL = min([cellfun('size', {trainD.data}, 2) cellfun('size', {testD.data}, 2)]);
                        for tr_ix = 1:length(trainD)
                            trainD(tr_ix).data = trainD(tr_ix).data(:, 1:dataL);
                        end; clear tr_ix
                        
                        for tr_ix = 1:length(testD)
                            testD(tr_ix).data = testD(tr_ix).data(:, 1:dataL);
                        end; clear tr_ix dataL
                        
                        thisTrainX = permute(cat(3, trainD.data), [3 1 2]);
                        thisTestX = permute(cat(3, testD.data), [3 1 2]);
                        
                        % Save to output
                        trainX{cv_ix} = thisTrainX;
                        testX{cv_ix} = thisTestX;
                        clear trainD testD
                    end; clear cv_ix
                end %if ~isTraj
                subTrainX{sfs_ix} = trainX;
                subTestX{sfs_ix} = testX;
            end
            
            %%
            subTrainX = cat(3, subTrainX{:});
            subTestX = cat(3, subTestX{:});
            %Merge sub feature sets back.
            trainX = cell(analysisParams.kcross, 1);
            testX = cell(analysisParams.kcross, 1);
            for cv_ix = 1:analysisParams.kcross
                trainX{cv_ix} = cat(3, subTrainX{cv_ix, :});
                testX{cv_ix} = cat(3, subTestX{cv_ix, :});
            end
            clear sfs_ix subTrainX subTestX subfsName cv_ix
            
            %% Classification
            for ml_ix = 1:length(mlNames)
                mlName = mlNames{ml_ix};
                predY = nan(size(trClass(eitherBool)));
                
                
                if strcmpi(mlName, 'rLDA')
                    addpath(bcipath{:});
                end
                
                %cross-validated
                for cv_ix = 1:analysisParams.kcross
                    
                    
                    nTrain = sum(trainBool(:, cv_ix));
                    nTest = sum(testBool(:, cv_ix));
                    
                    thisTrainY = trClass(trainBool(:, cv_ix));
                    thisTrainX = trainX{cv_ix};
                    thisTestX = testX{cv_ix};
                    [~, nNeurons, nTimes] = size(thisTrainX);
                    
                    % When working with neural trajectories...
                    %
                    % So far I don't have any ML methods that make explicit
                    % use of the temporal information, though I hope to achieve
                    % that with GLM (time factored out), or DBN/CRF.
                    %
                    % For the remaining methods (SVM, LDA, nbayes), we have to
                    % decide how to handle time.
                    if (isTraj || strcmpi(fsName, 'baseline,cue,delay,response'))...
                            && ~any(strcmpi({'GLM', 'DBN', 'CRF'}, mlName))
                        if false%any(strcmpi({'LDA'}, mlName))
                            %Option 1: Each time point is a different
                            %observation.
                            thisTrainX = reshape(permute(thisTrainX, [3 1 2]), [nTimes*nTrain nNeurons]);
                            thisTrainY = reshape(repmat(thisTrainY, 1, nTimes)', [], 1);
                            thisTestX = mean(thisTestX, 3);
                        else
                            %Option 2: Each time point is a different variable.
                            thisTrainX = reshape(thisTrainX, [nTrain, nNeurons*nTimes]);
                            thisTestX = reshape(thisTestX, [nTest, nNeurons*nTimes]);
                            %No adjustement to thisTrainY needed.
                        end
                    end
                    
                    % Do the machine learning.
                    if ~isempty(thisTrainX) && ~isempty(thisTestX) && ~isempty(thisTrainY)
                        
                        if doStandardizeX
                            x_bar = mean([thisTrainX;thisTestX]);
                            x_std = std([thisTrainX;thisTestX]);
                            
                            keepFeat = x_std > 0 & var(thisTrainX) > 0;
                            
                            thisTrainX = bsxfun(@rdivide, bsxfun(@minus, thisTrainX, x_bar),x_std);
                            thisTestX = bsxfun(@rdivide, bsxfun(@minus, thisTestX, x_bar), x_std);
                            
                            
                            thisTrainX = thisTrainX(:, keepFeat);
                            thisTestX = thisTestX(:, keepFeat);
                        end
                        
                        if strcmpi(mlName, 'SVM-Lin')
                            %Can I be clever about using time somehow?
                            %http://www.mathworks.com/help/stats/support-vector-machines-svm.html#bsr5b42
                            %Maybe libSVM is more straightforward?
                            t = templateSVM();
                            % 'svm' 'discriminant' 'naivebayes' 'knn' 'tree'
                            
                            model = fitcecoc(thisTrainX, thisTrainY,...
                                'Coding', 'onevsall', ...%'onevsone',...
                                'FitPosterior', true,...
                                'Learners', t,...
                                'CrossVal', 'off',...
                                'Prior', 'empirical',...
                                'Verbose', 0);
                            %kernel 'linearl', 'rbf', 'myCustomSigmoid'
                            predY(testBool(:, cv_ix)) = predict(model, thisTestX);
                            
                        elseif strcmpi(mlName, 'SVM-RBF')
                            t = templateSVM('KernelFunction', 'rbf');
                            model = fitcecoc(thisTrainX, thisTrainY,...
                                'Coding', 'onevsall', ...%'onevsone',...
                                'FitPosterior', true,...
                                'Learners', t,...
                                'CrossVal', 'off',...
                                'Prior', 'empirical',...
                                'Verbose', 0);
                            %kernel 'linearl', 'rbf', 'myCustomSigmoid'
                            predY(testBool(:, cv_ix)) = predict(model, thisTestX);
                            
                        elseif strcmpi(mlName, 'LDA')
                            %http://www.mathworks.com/help/stats/discriminant-analysis.html#btaf5dv
                            linclass = fitcdiscr(thisTrainX, thisTrainY);%, 'discrimType', 'quadratic');
                            predY(testBool(:, cv_ix)) = predict(linclass, thisTestX);
                            
                        elseif strcmpi(mlName, 'rLDA')
                            %BCILAB method.
                            model = ml_train({[ones(size(thisTrainX, 1), 1) thisTrainX], thisTrainY},...
                                {'lda',...
                                'lambda', [],...
                                'regularization', 'auto'});  %'auto' 'shrinkage' 'independence'
                            predictions = ml_predict([ones(size(thisTestX, 1), 1) thisTestX], model);
                            [~, predY(testBool(:, cv_ix))] = max(predictions{2}, [], 2);
                            
                        elseif strcmpi(mlName, 'NB')  % Naive Bayes
                            %My way: If we were to use pure spike counts then we
                            %would do 'model', 'poiss'.
                            %predY(~trainBool) = classifyNaiveBayesMulti(trainX, trainY, testX, 'model', 'norm');
                            model = fitcnb(thisTrainX, thisTrainY,...
                                'DistributionNames', 'kernel'); %,...% 'mn' might work for spike counts.
                            predY(testBool(:, cv_ix)) = predict(model, thisTestX);
                            
                        elseif strcmpi(mlName, 'GLM')
                            % Can I do a 2-way RM MANOVA? (i.e. 3-way RM ANOVA; random effects, nested within trial?)
                            % or GLM?  fRate = [class neuron time]; fitglm.
                            % Then...?
                            predY(testBool(:, cv_ix)) = nan;
                            
                        elseif strcmpi(mlName, 'DBN')
                            predY(testBool(:, cv_ix)) = nan;
                            %http://stats.stackexchange.com/questions/104697/how-do-i-detect-state-change-in-multivariate-time-series
                            %http://stackoverflow.com/questions/17487356/hidden-markov-model-for-multiple-observed-variables
                            %http://bnt.googlecode.com/svn/trunk/docs/usage_dbn.html
                            %Maybe DBN, maybe CRF
                            %Maybe all I have to do is take HMM's EM algorithm and tweak to handle multivariate emissions?
                            
                        end
                    end
                    clear nTrain nTest thisTrainY thisTrainX thisTestX model
                end; clear cv_ix
                
                predictedY{sess_ix, fs_ix, ml_ix}{pair_ix} = predY;
                
                if strcmpi(mlName, 'rLDA')
                    rmpath(bcipath{:});
                end
            end; clear ml_ix
        end; clear fs_ix
        
    end; clear pair_ix
    fprintf(['Finished ' num2str(sess_ix) ' of ' num2str(length(sessions)) '.\n']);
end
clear testBool ptb cDat trialBool Ynum Ystr nTrials this_sess

save(fullfile(paths.results, 'context_ml.mat'), ...
    'predictedY', 'actualY', 'actualYstr', 'analysisParams', 'dhParams',...
    'featureNames', 'mlNames', 'sessions', 'pairBoolOut');