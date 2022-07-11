% analyze_block_ml
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

sessions = my_sessions('Complementary');  %Need to convert more sessions for 'Blocks'
nSessions = length(sessions);
analysisParams = my_anaparams('classification');
%temporarily remove trajectory, to speed things up.
%%analysisParams.anaWins(strcmpi({analysisParams.anaWins.name}, 'trajectory')) = [];
nAna = length(analysisParams.anaWins);  % Number of analysis windows.

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
mlNames = {'SVM'};  %'SVM', 'LDA', 'NB', 'DNB'
nML = length(mlNames);

maxTrajDims = 8;
nTargPairs = 4;

%% Prepare output variables
accuracy_out = nan(nSessions, nFeatureSets, nML);

sess_ix = 1;
%%
for sess_ix = 1:nSessions
    %%
    this_sess = sessions(sess_ix);
    
    %% Load and process ptbmat
    ptb = load(fullfile(paths.preprocessed, 'ptb', this_sess.ptbfname), '-mat');
    [ptb.eData.trial, ptb.params] = getTrialStimInfo(ptb.eData.trial, ptb.params);
    ptb.eData.trial = getTrialBehavResult(ptb.eData.trial, ptb.params, analysisParams);
    ptb.eData.trial = getNewClass(ptb.eData.trial, analysisParams);
    
%     %% Calculate performance measures
%     %         cnfMat = [behav.TrueA behav.FalseA; ...
%     %                   behav.FalseB behav.TrueB];
%     %http://www.birmingham.ac.uk/Documents/college-les/psych/vision-laboratory/sdtintro.pdf
%     pmarg = 0.01; %If probability approaches 0 or 1, fix by this amount.
%     behav = procBehavBlocks(ptb.eData.trial);
%     tpr = behav.TrueA./(behav.TrueA + behav.FalseB);
%     tpr(tpr<pmarg)=pmarg; tpr(tpr>1-pmarg)=1-pmarg;
%     fpr = behav.FalseA./(behav.FalseA+behav.TrueB);
%     fpr(fpr<pmarg)=pmarg; fpr(fpr>1-pmarg)=1-pmarg;
%     ztpr = norminv(tpr, 0, 1);
%     zfpr = norminv(fpr, 0, 1);
%     d = ztpr - zfpr;
%     %     c = -(ztpr+zfpr)/2;
%     clear pmarg behav tpr fpr ztpr zfpr
    d = Inf*ones(size(ptb.eData.trial'));
    
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
    cDat.trial = cDat.trial(triageBool);
    d = d(triageBool);
    
    % Eliminate units that fire < 1Hz
    % Adds cDat.invalidUnits.
    % If thisParams.useUnits == 'merged' then it merges all units per
    % channel into one unit per channel.
    %53 units including noise.
    %19 sorted units after triage.
    [cDat, ~] = triageUnits(cDat, analysisParams);
    nUnits = size(cDat.trial(1).raster, 2);
    
    %% Dimensionality reduction on neural activity. Be agressive.
    
    % Done separately for each feature type
    
    %For now I will do feature reduction once on all trials, and not worry
    %about getting feature reduction parameters from training-only to apply
    %later to test data (especially important for dPCA and GPFA).
    
    fs_ix = 1; %find(strcmpi(featureNames, 'cue'));
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
        
        % For quick debug
        class_ix = 1;
        ml_ix = 1;
        cv_ix = 1;
        
        %% Classification
        allTrueY = [cDat.trial.newClass]';
        allPredY = nan(length(allTrueY), nML);
        
        %%
        thisBool = d'>1;
        thisX = X(thisBool, :, :);
        thisY = allTrueY(thisBool);
        testBool = cvDistrib(thisY, analysisParams.kcross, false);
        
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
                    model = fitcecoc(thisTrainX, thisTrainY,...
                        'Coding', 'onevsall', ...%'onevsone',...
                        'FitPosterior', true,...
                        'Learners', 'svm',...  % 'svm' 'discriminant' 'naivebayes' 'knn' 'tree'
                        'CrossVal', 'off',...
                        'Prior', 'empirical',...
                        'Verbose', 0);
                    %kernel 'linearl', 'rbf', 'myCustomSigmoid'
                    thisPredY = predict(model, thisTestX);
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
            allPredY(thisBool, ml_ix) = predY;
            
            clear predY
            if strcmpi(mlName, 'LDA')
                rmpath(bcipath{:});
            end
        end; clear ml_ix
        
        
        clear fsName trainX testX anaBool D
        match = bsxfun(@minus, allPredY, allTrueY) == 0;
        was_tested = ~any(isnan(allPredY), 2);
        if any(was_tested)
            accuracy_out(sess_ix, fs_ix, :) = sum(match(was_tested))./sum(was_tested);
        end
    end; clear fs_ix
    fprintf(['Finished ' num2str(sess_ix) ' of ' num2str(length(sessions)) '.\n']);
    
end
clear testBool ptb cDat trialBool Ynum Ystr nTrials this_sess

featureNames
mlNames
accuracy_out