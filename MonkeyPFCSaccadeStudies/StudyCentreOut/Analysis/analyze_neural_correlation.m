addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths DEF
addpath(genpath(fullfile(paths.ml3rd, 'DataHigh')));
addpath(fullfile(paths.ml3rd, 'dPCA', 'matlab'));
% addpath(genpath(fullfile(paths.ml3rd, 'L1General')));
% addpath(fullfile(paths.mlDevel, 'Helper'));
% addpath(fullfile(paths.ml3rd, 'dynamic_time_warping_v2.0'));
% addpath(fullfile(paths.ml3rd, 'DBA'));

%% Constants
sessions = my_sessions('CentreOut');
analysisParams = my_anaparams('classification');
analysisParams.classifyTarg = 'sacPol';
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

maxTrajDims = 8;
doStandardizeX = true;

%% Prepare output variables
predictedYnum = cell(length(sessions), nFeatureSets);
actualYnum = cell(length(sessions), 1);
rstat = cell(length(sessions), nFeatureSets);
targetClass = cell(length(sessions), 1);

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
    
%     Ynum = cat(1, cDat.trial.sacPol);  %tmp = Yang(:, 1)
%     Ynum = [sin(Ynum(:, 1)) cos(Ynum(:, 1))];
%     actualYnum{sess_ix} = Ynum;
    
    Ynum = cat(1, cDat.trial.sacEndXY);
    screenHalfsize = [512 384];%on 1024 x 768 screen
    Ynum = bsxfun(@rdivide, bsxfun(@minus, Ynum, screenHalfsize), screenHalfsize);
    actualYnum{sess_ix} = Ynum;
    
    targetClass{sess_ix} = [cDat.trial.newClass]';
    
    %% Prepare cross-validations so classes are well distributed across folds.
    % We will use these same CVs for feature reduction and for machine
    % learning.
    Ytemp = cat(1, cDat.trial.newClass);
    testBool = [cvDistrib(Ytemp, analysisParams.kcross) false(length(Ytemp), 1)];
    % Add one test cv which is really all trials are training.
    
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
                %newD(:, var(newD)==0) = [];
                
                % Put data into cell arrays for train and test.
                % Kind of ugly and redundant formatting, but necessary to match
                % format used by neural trajectories.
                
                for cv_ix = 1:size(testBool, 2)
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
                for cv_ix = 1:size(testBool, 2)
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
                        %convert to frate
                        nGoodSamps = sum(~isnan(raster), 3);
                        nSpikes = sum(raster > 0, 3);
                        frate = sqrt(1000*nSpikes./nGoodSamps);
                        thisTrainY = targetClass{sess_ix}(trainBool);
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
                        
                        if isempty(testIds)
                            testD = struct('data', {}, 'y', {});
                        else
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
                        trainX{cv_ix} = thisTrainX;  % Save to output
                        
                        if ~isempty(testD)
                            thisTestX = permute(cat(3, testD.data), [3 1 2]);
                            testX{cv_ix} = thisTestX;
                        end
                        
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
        fullX = cat(3, subTrainX{analysisParams.kcross+1, :});
        clear sfs_ix subTrainX subTestX subfsName cv_ix
        
        %% Multivariate regression
        
        [nTrials, nNeurons, nTimes] = size(fullX);
        dat = reshape(fullX, [nTrials, nNeurons*nTimes]);
        
        keepFeat = reshape(var(fullX) > 0, [1, nNeurons*nTimes]);
        dat = dat(:, keepFeat);
        
        if doStandardizeX
            dat = zscore(dat);
        end
        
        % Same as regress one at a time... which we need for r
        rstat{sess_ix, fs_ix} = nan(1, 2);
        lambda_out = nan(1, 2);
        for targ_dim = 1:2
            [B,FitInfo] = lasso(dat,Ynum(:, targ_dim), 'CV', 10);
            lambda_out(targ_dim) = FitInfo.LambdaMinMSE;
            pred_targ = dat*B(:, FitInfo.IndexMinMSE);
            r = corrcoef(Ynum(:, targ_dim), pred_targ);
            rstat{sess_ix, fs_ix}(targ_dim) = r(2,1);
        end
        
        % Cross-validated prediction
        predY = nan(size(Ynum));
        for cv_ix = 1:analysisParams.kcross
            trainBool = ~testBool(:, cv_ix);
            
            nTrain = sum(trainBool);
            nTest = sum(~trainBool);
            
            thisTrainY = Ynum(trainBool, :);
            thisTestY = Ynum(~trainBool, :);
            thisTrainX = trainX{cv_ix};
            thisTestX = testX{cv_ix};
            
            [~, nNeurons, nTimes] = size(thisTrainX);
            
            % So far I don't have a clever way of handling time as a factor
            if (isTraj || strcmpi(fsName, 'baseline,cue,delay,response'))
                if false
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
            
            if ~isempty(thisTrainX) && ~isempty(thisTestX) && ~isempty(thisTrainY)
                
                keepFeat = var(thisTrainX) > 0;
                if doStandardizeX
                    x_bar = mean([thisTrainX;thisTestX]);
                    x_std = std([thisTrainX;thisTestX]);
                    
                    keepFeat = keepFeat & x_std > 0;
                    
                    thisTrainX = bsxfun(@rdivide, bsxfun(@minus, thisTrainX, x_bar),x_std);
                    thisTestX = bsxfun(@rdivide, bsxfun(@minus, thisTestX, x_bar), x_std);
                end
                thisTrainX = thisTrainX(:, keepFeat);
                thisTestX = thisTestX(:, keepFeat);
                
                for targ_dim = 1:2
                    [B,FitInfo] = lasso(thisTrainX, thisTrainY(:, targ_dim),...
                        'Lambda', lambda_out(targ_dim));
                    predY(~trainBool, targ_dim) = thisTestX*B;
                end
            
            end
        end
        predictedYnum{sess_ix, fs_ix} = predY;
        clear predY
        clear fsName trainX testX anaBool D

    end; clear fs_ix
    
    fprintf(['Finished ' num2str(sess_ix) ' of ' num2str(length(sessions)) '.\n']);
    clear testBool ptb cDat trialBool Ynum Ystr nTrials this_sess
end
save(fullfile(paths.results, 'correlation.mat'), ...
    'actualYnum', 'analysisParams', 'dhParams',...
    'featureNames', 'predictedYnum', 'sessions', 'targetClass');