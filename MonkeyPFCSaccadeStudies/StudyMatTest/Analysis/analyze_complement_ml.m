% analyze_complement_ml
% A script to
%%
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths DEF

% %% BCILAB - startup is very slow.
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

%% Paths and Constants
addpath(genpath(fullfile(paths.ml3rd, 'DataHigh')));

sessions = my_sessions('Complementary');
nSessions = length(sessions);
analysisParams = my_anaparams('classification');
nAna = length(analysisParams.anaWins);  % Number of analysis windows.
analysisParams.anaWins(strcmpi({analysisParams.anaWins.name}, 'trajectory')) = [];

% Feature reduction parameters (using datahigh)
dhParams = struct(...
    'binWidth', analysisParams.binWidth,...
    'use_sqrt', true,...
    'kern', analysisParams.kernSD,...
    'trial_average', false,...
    'keep_neurons', true(1,1));

% Feature set names
frNames = {'smooth', 'gpfa'}; % Feature reduction names for anaWins.name == trajectory
featureNames = {analysisParams.anaWins.name};
traj_bool = strcmpi(featureNames, 'trajectory');
if any(traj_bool)
    featureNames(traj_bool) = [];
    featureNames = [featureNames frNames];
end
featureNames = [featureNames 'baseline,visual,cue,delay']; %Concatenated spike counts.
nFeatureSets = length(featureNames);
clear traj_bool

% Machine learning names.
%CCA has no advantage over LDA for spike counts.
%CCA can eliminate time as a factor when using trajectories.
mlNames = {'SVM'};  %, 'LDA', 'NB', 'DNB'
nML = length(mlNames);

maxTrajDims = 8;
nTargPairs = 4;

%% Prepare output variables
accuracy_out = cell(length(sessions), nFeatureSets, nML);

sess_ix = 1;
%%
for sess_ix = 1:nSessions
    %%
    this_sess = sessions(sess_ix);
    
    %% Load and process ptbmat
    ptb = load(fullfile(paths.preprocessed, 'ptb', this_sess.ptbfname), '-mat');
    [ptb.eData.trial, ptb.params] = getTrialStimInfo(ptb.eData.trial, ptb.params);
    ptb.eData.trial = getTrialBehavResult(ptb.eData.trial, ptb.params, analysisParams);
    pairBool = getComplementaryMembership(ptb.eData.trial);
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
    triageBool = triageTrials(cDat, analysisParams);
    pairBool = pairBool(triageBool, :, :, :);
    trialBoolOut{sess_ix} = triageBool;
    cDat.trial = cDat.trial(triageBool); clear trialBool
    
    % Eliminate units that fire < 1Hz
    % Adds cDat.invalidUnits.
    % If thisParams.useUnits == 'merged' then it merges all units per
    % channel into one unit per channel.
    %53 units including noise.
    %19 sorted units after triage.
    [cDat, fullInvalidOut{sess_ix}] = triageUnits(cDat, analysisParams);
    nUnits = size(cDat.trial(1).raster, 2);
    
    % The targets and colours we can test with this session
    [availTarg, availCol] = find(squeeze(~any(sum(pairBool)==0, 3)));
    
    %% Dimensionality reduction on neural activity. Be agressive.
    
    % Done separately for each feature type
    
    %For now I will do feature reduction once on all trials, and not worry
    %about getting feature reduction parameters from training-only to apply
    %later to test data (especially important for dPCA and GPFA).
    
    fs_ix = find(strcmpi(featureNames, 'cue'));
    %%
    for fs_ix = 1:nFeatureSets
        %%
        fsName = featureNames{fs_ix};
        isTraj = any(strcmpi(frNames, fsName));
        subfsName = strsplit(fsName, ',');
        subX = cell(1, length(subfsName));
        for sfs_ix = 1:length(subfsName)
            sfsName = subfsName{sfs_ix};
            
            % Convert to DataHigh format
            if isTraj
                ana_win = analysisParams.anaWins(strcmpi({analysisParams.anaWins.name}, 'trajectory'));
            else
                ana_win = analysisParams.anaWins(strcmpi({analysisParams.anaWins.name}, sfsName));
            end
            cDat = trimTrials(cDat, ana_win);
            D = cDat2DataHigh(cDat.trial);
            
            % Reset DataHigh parameters
            dhParams.binWidth = analysisParams.binWidth;
            dhParams.kern = analysisParams.kernSD;
            dhParams.keep_neurons = cDat.invalidUnits == DEF.unitState.valid;
            nDimRed = sum(dhParams.keep_neurons);
            dhRedMode = -1;
            
            if ~isTraj
                lcdLength = nan(size(D));
                for tr_ix = 1:length(D)
                    lcdLength(tr_ix) = size(D(tr_ix).data, 2);
                end; clear tr_ix
                dhParams.binWidth = min(lcdLength);
                dhParams.kern = 0;
            elseif strcmpi(sfsName, 'pca')
                dhRedMode = 1;
                nDimRed = maxTrajDims;
            elseif strcmpi(sfsName, 'gpfa')
                dhRedMode = 5;
                nDimRed = maxTrajDims;
            end
            [newD, ~, ~] = reducedims(D, dhRedMode, nDimRed, dhParams);
            
            % Maybe necessary to trim trials to same length.
            dataL = min(cellfun('size', {newD.data}, 2));
            for tr_ix = 1:length(newD)
                newD(tr_ix).data = newD(tr_ix).data(:, 1:dataL);
            end; clear tr_ix
            
            timeVec = cDat.trial(1).tVec(cDat.trial(1).timeBool);
            if mod(length(timeVec), dhParams.binWidth) > 0
                timeVec = [timeVec nan(1, dhParams.binWidth - mod(length(timeVec), dhParams.binWidth))];
            end
            timeVec = nanmean(reshape(timeVec, dhParams.binWidth, []));
            timeVec = timeVec(1:dataL);
            
            subX{1, sfs_ix} = permute(cat(3, newD.data), [3 1 2]);
        end
        X = cat(3, subX{:});
        
        %% TODO: Standardize X by baselineX, PER BLOCK!
        %For now, just divide by baseline
        block_id = [cDat.trial.cueTargRuleGroup];
        block_uq = unique(block_id);
        if strcmpi(fsName, 'baseline')
            baselineX = X;
            bl_avg = nan(length(block_uq), nUnits);
            bl_std = nan(length(block_uq), nUnits);
            for block_ix = 1:length(block_uq)
                block_bool = block_id == block_uq(block_ix);
                bl_avg(block_ix, :) = squeeze(mean(baselineX(block_bool, :)));
                bl_std(block_ix, :) = squeeze(std(baselineX(block_bool, :)));
            end
        elseif any(strcmpi(fsName, {'gpfa', 'dpca'}))
            bl_timeEdges = analysisParams.anaWins(strcmpi({analysisParams.anaWins.name}, 'baseline')).winEdges;
            blTimeBool = timeVec >= bl_timeEdges(1) & timeVec <= bl_timeEdges(2);
            traj_bl_avg = nan(length(block_uq), maxTrajDims);
            traj_bl_std = nan(length(block_uq), maxTrajDims);
            for block_ix = 1:length(block_uq)
                block_bool = block_id == block_uq(block_ix);
                X_bl = mean(X(block_bool, :, blTimeBool), 3);
                traj_bl_avg(block_ix, :) = mean(X_bl);
                traj_bl_std(block_ix, :) = std(X_bl);
            end
        end
        
        for block_ix = 1:length(block_uq)
            block_bool = block_id == block_uq(block_ix);
            if ~isTraj
                X(block_bool, :, :) = bsxfun(@rdivide,...
                    X(block_bool, :, :), bl_avg(block_ix, :));
                    %bsxfun(@minus, X(block_bool, :, :), bl_avg(block_ix, :)),...
                    %                         bl_std(block_ix, :));
            else
                X(block_bool, :, :) = bsxfun(@rdivide,...
                    X(block_bool, :, :), traj_bl_avg(block_ix, :));
                    %bsxfun(@minus, X(block_bool, :, :), traj_bl_avg(block_ix, :)),...
                    %                         traj_bl_std(block_ix, :));
            end
        end
        
        % For quick debug
        class_ix = 1;
        ml_ix = 1;
        cv_ix = 1;
        
        %% Classification
        for class_ix = 1:length(availTarg)
            thisBool = squeeze(pairBool(:, availTarg(class_ix), :, availCol(class_ix)));
            thisX = X(any(thisBool, 2), :, :);
            thisY = -1*double(thisBool(:, 1)) + 1*double(thisBool(:, 2));
            thisY = thisY(any(thisBool, 2));
            
            %TODO: Option to match nTraining
            testBool = cvDistrib(thisY, analysisParams.kcross, true);
            
            
            %%
            for ml_ix = 1:nML
                %%
                mlName = mlNames{ml_ix};
                
                if strcmpi(mlName, 'LDA')
                    addpath(bcipath{:});
                end
                
                predY = nan(size(thisY));
                %% cross-validated
                for cv_ix = 1:analysisParams.kcross
                    %% Prepare the training and testing X & Y
                    trainBool = ~testBool(:, cv_ix);
                    nTrain = sum(trainBool);
                    nTest = sum(~trainBool);
                    thisTrainX = thisX(trainBool, :, :);
                    thisTrainY = thisY(trainBool);
                    thisTestX = thisX(~trainBool, :, :);
                    thisTestY = thisY(~trainBool);
                    [~, nNeurons, nTimes] = size(thisTrainX);
                    
                    % Handle the 'time' dimension.
                    if (isTraj || length(subfsName) > 1)...
                            && any(strcmpi({'SVM', 'LDA', 'NB'}, mlName))
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
                    
                    %% Do the machine learning.
                    if strcmpi(mlName, 'SVM')
                        %Can I be clever about using time somehow?
                        %http://www.mathworks.com/help/stats/support-vector-machines-svm.html#bsr5b42
                        %Maybe libSVM is more straightforward?
                        model = fitcsvm(thisTrainX,thisTrainY,...
                            'Standardize', 'on',...
                            'KernelFunction', 'linear',...
                            'KernelScale', 'Auto');
                        [thisPredY, scores] = predict(model, thisTestX);
                    elseif strcmpi(mlName, 'LDA')
                        %                     %http://www.mathworks.com/help/stats/discriminant-analysis.html#btaf5dv
                        %                     linclass = fitcdiscr(thisTrainX, thisTrainY);%, 'discrimType', 'quadratic');
                        %                     predY(~trainBool) = predict(linclass, thisTestX);
                        
                        %BCILAB method.
                        model = ml_train({[ones(size(thisTrainX, 1), 1) thisTrainX], thisTrainY},...
                            {'lda',...
                            'lambda', [],...
                            'regularization', 'auto'});  %'auto' 'shrinkage' 'independence'
                        predictions = ml_predict([ones(size(thisTestX, 1), 1) thisTestX], model);
                        [~, thisPredY] = max(predictions{2}, [], 2);
                    elseif strcmpi(mlName, 'NB')  % Naive Bayes
                        %My way: If we were to use pure spike counts then we
                        %would do 'model', 'poiss'.
                        %predY(~trainBool) = classifyNaiveBayesMulti(trainX, trainY, testX, 'model', 'norm');
                        model = fitcnb(thisTrainX, thisTrainY,...
                            'DistributionNames', 'kernel'); %,...% 'mn' might work for spike counts.
                        thisPredY = predict(model, thisTestX);
                    elseif strcmpi(mlName, 'GLM')
                        % Can I do a 2-way RM MANOVA? (i.e. 3-way RM ANOVA; random effects, nested within trial?)
                        % or GLM?  fRate = [class neuron time]; fitglm.
                        % Then...?
                        thisPredY = nan(sum(~trainBool), 1);
                    elseif strcmpi(mlName, 'DNB')
                        thisPredY = nan(sum(~trainBool), 1);
                        %http://stats.stackexchange.com/questions/104697/how-do-i-detect-state-change-in-multivariate-time-series
                        %http://stackoverflow.com/questions/17487356/hidden-markov-model-for-multiple-observed-variables
                        %http://bnt.googlecode.com/svn/trunk/docs/usage_dbn.html
                        %Maybe DNB, maybe CRF
                        %Maybe all I have to do is take HMM's EM algorithm and tweak to handle multivariate emissions?
                    end
                    predY(~trainBool) = thisPredY;
                    clear nTrain nTest thisTrainY thisTrainX thisTestX model trainBool
                end; clear cv_ix
                %%
                accuracy_out{sess_ix, fs_ix, ml_ix}(class_ix) = sum(predY == thisY)/length(predY);
                clear predY
                if strcmpi(mlName, 'LDA')
                    rmpath(bcipath{:});
                end
            end; clear ml_ix
            
        end; clear class_ix
        clear fsName trainX testX anaBool D
    end; clear fs_ix
    fprintf(['Finished ' num2str(sess_ix) ' of ' num2str(length(sessions)) '.\n']);
    featureNames
    cat(1, accuracy_out{sess_ix, :, 1})'
end
clear testBool ptb cDat trialBool Ynum Ystr nTrials this_sess

featureNames
for sess_ix = 1:size(accuracy_out, 1)
    cat(1, accuracy_out{sess_ix, :, 1})'
end




