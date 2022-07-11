addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths DEF
addpath(genpath(fullfile(paths.ml3rd, 'DataHigh')));
addpath(fullfile(paths.ml3rd, 'dPCA', 'matlab'));
% addpath(genpath(fullfile(paths.ml3rd, 'L1General')));
% addpath(fullfile(paths.mlDevel, 'Helper'));
% addpath(genpath(fullfile(paths.ml3rd, 'CCToolbox')));
% addpath(fullfile(paths.ml3rd, 'dynamic_time_warping_v2.0'));
% addpath(fullfile(paths.ml3rd, 'DBA'));
% addpath(genpath(fullfile(paths.mlDevel, 'MachineLearning')));
% addpath(genpath(fullfile(paths.mlDevel, 'Evaluate')));
%%
% %% BCILAB - startup is very slow.
% %TODO: Get rid of this once I find better rLDA and rLogReg.
prevwd = pwd();
bcipath = fullfile(paths.ml3rd,'BCILAB');
addpath(bcipath);
bcilab;
cd(prevwd);
p = strsplit(path, pathsep);
isbcipath = false(size(p));
for p_ix = 1:length(p)
    isbcipath(p_ix) = length(p{p_ix}) >= length(bcipath) && strcmpi(p{p_ix}(1:length(bcipath)), bcipath);
end
bcipath = p(isbcipath);
rmpath(bcipath{:});
clear prevwd p isbcipath p_ix CURRENTSTUDY

%% Constants
sessions = my_sessions('CentreOut');
analysisParams = my_anaparams('classification');
nAna = length(analysisParams.anaWins);

tuningParams = my_anaparams('tuning');  % Will be used for zerovar epochs

% Additional analysis window structure needed for canonical trajectories
epochAnaWin.cue = analysisParams.anaWins(strcmpi({analysisParams.anaWins.name}, 'cue'));
epochAnaWin.delay = analysisParams.anaWins(strcmpi({analysisParams.anaWins.name}, 'delay'));
epochAnaWin.resp = analysisParams.anaWins(strcmpi({analysisParams.anaWins.name}, 'response'));

% Feature reduction parameters (using datahigh)
dhParams = struct(...
    'binWidth', analysisParams.binWidth,...
    'use_sqrt', true,...
    'kern', analysisParams.kernSD,...
    'trial_average', false,...
    'keep_neurons', true(1,1));

% Feature set names
featureNames = {analysisParams.anaWins.name};
% baseline cue delay response fullTrial trajectory

% Feature reduction names for anaWins.name == trajectory
trajNames = {'smooth' 'pca', 'fa', 'gpfa', 'dPCA', 'canon_cue', 'canon_delay', 'canon_resp'};
traj_bool = strcmpi(featureNames, 'trajectory');
if any(traj_bool)
    featureNames(traj_bool) = [];
    featureNames = [featureNames trajNames];
end
featureNames = [featureNames 'baseline,cue,delay,response']; %Concatenated spike counts.
clear traj_bool

nFeatureSets = length(featureNames);

% Machine learning names.
mlNames = {'SVM-Lin', 'SVM-RBF', 'LDA', 'rLDA'};  %DBN, DTW, NB
nML = length(mlNames);

maxTrajDims = 8;
doStandardizeX = true;

%% Prepare output variables
predictedY = cell(length(sessions), nFeatureSets, nML);
actualY = cell(length(sessions), 1);
actualYstr = actualY;

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
    
    % Eliminate units that...
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
    % We will use these same CVs for feature reduction and for machine
    % learning.
    testBool = cvDistrib(Ynum, analysisParams.kcross);
    
    % baseline cue delay response smooth pca gpfa dPCA canon_cue canon_delay canon_resp baseline,cue,delay
    %fs_ix = find(strcmpi(featureNames, 'canon_resp'));
    
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
    
    %% For each Feature Set
    for fs_ix = 1:nFeatureSets
        fsName = featureNames{fs_ix};
        isTraj = any(strcmpi(trajNames, fsName));
        subfsNames = strsplit(fsName, ',');
        
        subTrainX = cell(1, length(subfsNames));
        subTestX = cell(1, length(subfsNames));
        
        %%
        for sfs_ix = 1:length(subfsNames)
            %% Do feature reduction
            sfsNames = strsplit(subfsNames{sfs_ix}, '_'); %split again in case canonical
            sfsName = sfsNames{1};
            
            % Convert the data for dataHigh
            if isTraj
                ana_win = analysisParams.anaWins(strcmpi({analysisParams.anaWins.name}, 'trajectory'));
            else
                ana_win = analysisParams.anaWins(strcmpi({analysisParams.anaWins.name}, sfsName));
            end
            cDat = trimTrials(cDat, ana_win);
            D = cDat2DataHigh(cDat.trial);
            clear ana_win
            
            %Remove neurons with one or more bad epochs
            for tr_ix = 1:length(D)
                D(tr_ix).data = D(tr_ix).data(keepBool, :);
            end
            
            % Reset the DataHigh parameters
            dhParams.binWidth = analysisParams.binWidth;
            dhParams.kern = analysisParams.kernSD;
            dhParams.keep_neurons = cDat.invalidUnits(keepBool) == DEF.unitState.valid;
            nDimRed = sum(dhParams.keep_neurons);
            dhRedMode = -1;
            
            trainX = cell(analysisParams.kcross, 1);
            testX = trainX;
            
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
                
                % Trimming trials (esp response epoch) may result in 0-var
                % for some units.
                % Not necessary to clean here. It'll happen before ML.
                % newD(:, var(newD)==0) = [];
                
                % Put data into cell arrays for train and test.
                % Kind of ugly and redundant formatting, but necessary to match
                % format used by neural trajectories.
                
                for cv_ix = 1:analysisParams.kcross
                    trainBool = ~testBool(:, cv_ix);
                    trainX{cv_ix} = newD(trainBool, :, :);
                    testX{cv_ix} = newD(~trainBool, :, :);
                end; clear cv_ix
                
                clear dhRedMode lcdLength nDimRed newD dataL trainBool
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
                elseif strcmpi(sfsName, 'dPCA')
                    %keep -1 for smooth. TODO: Make dPCA part of DataHigh.
                elseif strcmpi(sfsName, 'canon')
                    %keep -1 for smooth.
                end
                
                % Trajectory transformations must be calculated using training
                % data only.
                for cv_ix = 1:analysisParams.kcross
                    abort_cv = false;
                    trainBool = ~testBool(:, cv_ix);
                    
                    %% Get the trajectories using training data only.
                    [trainD, C, lat, dhRedParams] = reducedims(D(trainBool), dhRedMode, nDimRed, dhParams);
                    
                    % Additional steps required for dPCA, canon
                    if strcmpi(fsName, 'dPCA')
                        [trainD, W, V, whichMarg, dPCAparams] = dPCAofDataHigh(trainD,...
                            'margNames', {'Location-Time', 'Time'}, ...
                            'margCombs', {{1, [1 2]}, {2}});
                        %optimLambda
                        %Find up to maxTrajDims best Location-Time components
                        marg_ix = find(strcmpi(dPCAparams.margNames, 'Location-Time'));
                        dOut = min(maxTrajDims, sum(whichMarg == marg_ix));
                        comp_ids = find(whichMarg == marg_ix, dOut, 'first');
                        W = W(:, comp_ids);
                        trainD = centerAndTransformStructArray(trainD, W);
                        clear V whichMarg dPCAparams margNames margCombs
                        clear marg_ix dOut comp_ids
                    elseif strcmpi(sfsName, 'canon')
                        raster = getEpochRaster(cDat, epochAnaWin.(sfsNames{2}));
                        raster = raster(trainBool, keepBool, :);
                        %Trim off times with nan so all trials may have the same number
                        %of timepoints
                        %convert to frate
                        nGoodSamps = sum(~isnan(raster), 3);
                        nSpikes = sum(raster > 0, 3);
                        frate = sqrt(1000*nSpikes./nGoodSamps);
                        thisTrainY = Ynum(trainBool);
                        try
                            [~, ~, stats] = manova1(frate, thisTrainY);
                            W = stats.eigenvec(:, 1:maxTrajDims);
                            trainD = centerAndTransformStructArray(trainD, W);
                        catch ME
                            abort_cv = true;
                            if strcmpi(ME.identifier, 'stats:manova1:SingularSumSquares2')
                                fprintf('MANOVA1 failed for session %i %s\n', sess_ix, sfsNames{2});
                            else
                                rethrow(ME)
                            end
                        end
                        clear raster goodBool spkCount thisTrainY stats
                    end
                    
                    %% Get the test trajectories using the training loading mat C.
                    if ~abort_cv
                        % Trial info
                        thisBool = ~trainBool;  %test trials
                        testIds = find(thisBool);
                        tId = num2cell(testIds);
                        
                        %Convert testing dat to transformable format
                        %using reducedims with dhRedMode == -1
                        [testD, ~, ~, ~] = reducedims(D(thisBool), -1, size(D(1).data, 1), dhParams);
                        [testD.y] = deal(testD.data);  % Copy data into y.
                        t_length = nan(size(testD));
                        for tr_ix = 1:numel(testD)
                            t_length(tr_ix) = size(testD(tr_ix).data, 2);
                        end
                        t_length = num2cell(t_length);
                        [testD.trialId] = tId{:};
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
                        elseif any(strcmpi({'dPCA', 'canon'}, sfsName))
                            testD = centerAndTransformStructArray(testD, W);
                        end
                        
                        %diag(corr([trainD.data]', [testD.data]'))
                        clear trainBool thisBool testIds tId s C lat d
                        
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
                    end
                end; clear cv_ix  % for cv_ix
            end
            
            subTrainX{sfs_ix} = trainX;
            subTestX{sfs_ix} = testX;
            clear dhRedMode nDimRed dhRedParams trainX testX
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
            %%
            % Put the ml_ix on the outer loop in case we have to change any
            % paths.
            
            mlName = mlNames{ml_ix};
            predY = nan(size(Ynum));
            
            if strcmpi(mlName, 'rLDA')
                addpath(bcipath{:});
            end
            
            %cross-validated
            for cv_ix = 1:analysisParams.kcross
                trainBool = ~testBool(:, cv_ix);
                
                nTrain = sum(trainBool);
                nTest = sum(~trainBool);
                
                thisTrainY = Ynum(trainBool);
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
                        predY(~trainBool) = predict(model, thisTestX);
                        
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
                        predY(~trainBool) = predict(model, thisTestX);
                        
                    elseif strcmpi(mlName, 'LDA')
                        %http://www.mathworks.com/help/stats/discriminant-analysis.html#btaf5dv
                        linclass = fitcdiscr(thisTrainX, thisTrainY);%, 'discrimType', 'quadratic');
                        predY(~trainBool) = predict(linclass, thisTestX);
                        
                    elseif strcmpi(mlName, 'rLDA')
                        %BCILAB method.
                        model = ml_train({[ones(size(thisTrainX, 1), 1) thisTrainX], thisTrainY},...
                            {'lda',...
                            'lambda', [],...
                            'regularization', 'auto'});  %'auto' 'shrinkage' 'independence'
                        predictions = ml_predict([ones(size(thisTestX, 1), 1) thisTestX], model);
                        [~, predY(~trainBool)] = max(predictions{2}, [], 2);
                        
                    elseif strcmpi(mlName, 'NB')  % Naive Bayes
                        %My way: If we were to use pure spike counts then we
                        %would do 'model', 'poiss'.
                        %predY(~trainBool) = classifyNaiveBayesMulti(trainX, trainY, testX, 'model', 'norm');
                        model = fitcnb(thisTrainX, thisTrainY,...
                            'DistributionNames', 'kernel'); %,...% 'mn' might work for spike counts.
                        predY(~trainBool) = predict(model, thisTestX);
                        
                    elseif strcmpi(mlName, 'GLM')
                        % Can I do a 2-way RM MANOVA? (i.e. 3-way RM ANOVA; random effects, nested within trial?)
                        % or GLM?  fRate = [class neuron time]; fitglm.
                        % Then...?
                        predY(~trainBool) = nan;
                        
                    elseif strcmpi(mlName, 'DBN')
                        predY(~trainBool) = nan;
                        %http://stats.stackexchange.com/questions/104697/how-do-i-detect-state-change-in-multivariate-time-series
                        %http://stackoverflow.com/questions/17487356/hidden-markov-model-for-multiple-observed-variables
                        %http://bnt.googlecode.com/svn/trunk/docs/usage_dbn.html
                        %Maybe DBN, maybe CRF
                        %Maybe all I have to do is take HMM's EM algorithm and tweak to handle multivariate emissions?
                        
                    end
                end
                clear nTrain nTest thisTrainY thisTrainX thisTestX model trainBool
            end; clear cv_ix
            
            predictedY{sess_ix, fs_ix, ml_ix} = predY;
            clear predY
            
            if strcmpi(mlName, 'rLDA')
                rmpath(bcipath{:});
            end
            
        end; clear ml_ix
        clear fsName trainX testX anaBool D
        
    end; clear fs_ix
    
    fprintf(['Finished ' num2str(sess_ix) ' of ' num2str(length(sessions)) '.\n']);
    clear testBool ptb cDat trialBool Ynum Ystr nTrials this_sess
end
save(fullfile(paths.results, 'classification.mat'), ...
    'actualY', 'actualYstr', 'analysisParams', 'dhParams',...
    'featureNames', 'mlNames', 'predictedY', 'sessions');