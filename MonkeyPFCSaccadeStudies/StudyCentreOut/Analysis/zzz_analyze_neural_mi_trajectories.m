% This script analyzes the mutual information among multi-dimensional
% neural trajectory and target location at each time-step in the trial,
% and compares it to a null-MI (1000 shuffles).

%% Paths and Constants
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths DEF

addpath(fullfile(paths.ml3rd, 'MIToolbox'));
addpath(genpath(fullfile(paths.ml3rd, 'DataHigh')));
addpath(fullfile(paths.ml3rd, 'dPCA', 'matlab'));

sessions = my_sessions('CentreOut');
nSessions = length(sessions);

analysisParams = my_anaparams('trajOnly');
nAna = length(analysisParams.anaWins);

epochAnaParams = my_anaparams('tuning');
baseAnaWin = epochAnaParams.anaWins(strcmpi({epochAnaParams.anaWins.name}, 'baseline'));
delayAnaWin = epochAnaParams.anaWins(strcmpi({epochAnaParams.anaWins.name}, 'delay'));

% Feature reduction parameters (using datahigh)
dhParams = struct(...
    'binWidth', analysisParams.binWidth,...
    'use_sqrt', true,...
    'kern', analysisParams.kernSD,...
    'trial_average', false,...
    'keep_neurons', true(1,1));
featureNames = {'smooth' 'pca', 'fa', 'gpfa', 'dPCA', 'canon_delay'};
nFeatureSets = length(featureNames);
maxTrajDims = 8;

trajMiStatP = cell(nSessions, nAna, nFeatureSets);%Preallocate
baseMiStatP = nan(nSessions, 2);

%%
for sess_ix = 1:nSessions
    %%
    this_sess = sessions(sess_ix);
    thisParams = analysisParams;
    
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
    cDat.trial = cDat.trial(trialBool); clear trialBool
    
    % Eliminate units that fire < 1Hz
    % Adds cDat.invalidUnits.
    % If thisParams.useUnits == 'merged' then it merges all units per
    % channel into one unit per channel.
    cDat = triageUnits(cDat, analysisParams);
    
    %% Easy-access, map classes -> directions
    classId = [cDat.trial.newClass]';
    [uqClasses, inClassIx] = unique(classId);
    nClasses = length(uqClasses);
    classStr = {cDat.trial(inClassIx).newClassStr};
    classDir = cat(1, cDat.trial(inClassIx).targPol); classDir = classDir(:, 1);
    clear inClassIx
    
    nTrials = length(cDat.trial);
    
    %% For each analysis window
    for win_ix = 1:nAna
        
        ana_win = analysisParams.anaWins(win_ix);
        
        %% For each feature set
        %fs_ix = find(strcmpi(featureNames, 'canon_delay'));
        for fs_ix = 1:nFeatureSets
            fsName = featureNames{fs_ix};
            
            %% Trim trials, depending on featureset
            ana_win.commonWin = strcmpi(fsName, 'dPCA');  % dPCA requires even-length trials
            cDat = trimTrials(cDat, ana_win);
            D = cDat2DataHigh(cDat.trial);
            
            if strcmpi(fsName, 'dPCA')
                tVec = cDat.trial(1).tVec(cDat.trial(1).timeBool);
            else
                tVec = ana_win.winEdges(1):ana_win.winEdges(2);
            end
            tVec = tVec(1:end-mod(length(tVec), analysisParams.binWidth));
            tVec = reshape(tVec, analysisParams.binWidth, []);
            step_edges = [tVec(1,:);tVec(end,:)];
            tVec = mean(tVec);
            nSteps = length(tVec);
            
            %% Map each trial's trajectory-bin to the common steps
            % This will be used to align data during MI calculation
            nSamples = cell2mat(cellfun(@size, {D.data},...
                repmat({2}, 1, length(D)), 'UniformOutput', false));
            if strcmpi(fsName, 'dPCA')
                %1:1 because previously trimmed to be common
                tr_first_steps = ones(1, nTrials);
            else
                tr_first_steps = nan(1, nTrials);
                for tr_ix = 1:nTrials
                    tr_first_steps(tr_ix) = find(tVec > cDat.trial(tr_ix).tVec(1), 1, 'first');
                end
            end
            tr_last_steps = tr_first_steps - 1 + floor(nSamples/analysisParams.binWidth);
                        
            %% Get the trajectories
            
            % Reset the DataHigh parameters
            dhParams.binWidth = analysisParams.binWidth;
            dhParams.kern = analysisParams.kernSD;
            dhParams.keep_neurons = cDat.invalidUnits == DEF.unitState.valid;
            nDimRed = sum(dhParams.keep_neurons);
            
            dhRedMode = -1;
            %smooth: -1; pca:1; fa:3; lda: 4, gpfa:5
            if strcmpi(fsName, 'pca')
                dhRedMode = 1;
                nDimRed = maxTrajDims;
            elseif strcmpi(fsName, 'fa')
                dhRedMode = 3;
                nDimRed = maxTrajDims;
            elseif strcmpi(fsName, 'lda')
                dhRedMode = 4;
                nDimRed = length(unique({D.condition}));
            elseif strcmpi(fsName, 'gpfa')
                dhRedMode = 5;
                dhParams.kern = 0;
                nDimRed = maxTrajDims;
            elseif strcmpi(fsName, 'dPCA')
                %keep -1 for smooth. TODO: Make dPCA part of DataHigh.
            elseif strcmpi(fsName, 'canon_delay')
                %keep -1 for smooth. TODO: Compare manova1 result with that
                %from datahigh for LDA when using only baseline epoch
            end
            
            [trajectories, C, lat, dhRedParams] = reducedims(D, dhRedMode, nDimRed, dhParams);
            
            if strcmpi(fsName, 'dPCA')
                [trajectories, W, V, whichMarg, dPCAparams] = dPCAofDataHigh(trajectories,...
                            'margNames', {'Location-Time', 'Time'}, ...
                            'margCombs', {{1, [1 2]}, {2}});
                %Find up to maxTrajDims best Location-Time components
                marg_ix = find(strcmpi(dPCAparams.margNames, 'Location-Time'));
                dOut = min(maxTrajDims, sum(whichMarg == marg_ix));
                comp_ids = find(whichMarg == marg_ix, dOut, 'first');
                W = W(:, comp_ids);
                trajectories = centerAndTransformStructArray(trajectories, W);
            elseif strcmpi(fsName, 'canon_delay')
                %Get delay epoch activity
                raster = getEpochRaster(cDat, delayAnaWin);
                %Trim off times with nan so all trials may have the same.
                goodBool = ~squeeze(any(any(isnan(raster), 2), 1));
                
                spkCount = sqrt(sum(raster(:, :, goodBool), 3));
                [~, ~, stats] = manova1(spkCount, classId);
                trajectories = centerAndTransformStructArray(trajectories,...
                    stats.eigenvec(:, 1:maxTrajDims));
            end
            nTrajDim = size(trajectories(1).data, 1);
            
            %% Calculate MI at each step
            trajMiStatP{sess_ix, win_ix, fs_ix} = nan(nSteps, 3);
            trajMiStatP{sess_ix, win_ix, fs_ix}(:,1) = tVec';
            for step_ix = 1:nSteps
                %Figure out which trials have data for this step.
                tr_bool = step_ix >= tr_first_steps & step_ix <= tr_last_steps;
                tr_inds = find(tr_bool);
                %Get the data using only these trials
                step_dat = nan(sum(tr_bool), nTrajDim);
                for tr_ix = 1:length(tr_inds)
                    step_dat(tr_ix, :) = trajectories(tr_inds(tr_ix)).data(:, tr_first_steps(tr_inds(tr_ix))-1+step_ix);
                end
                step_class = classId(tr_bool);
                [class_count, class_id] = hist(step_class, unique(step_class));
                
                if ~isempty(class_id) && ~any(class_count < analysisParams.minTrialsPerGroup)
                    %Calculate MI
                    thisMi = mi(step_dat, step_class);
                    %Calculate shuffled-label MI
                    miPerm = nan(1, analysisParams.statPerm);
                    for shuf_ix = 1:analysisParams.statPerm
                        miPerm(shuf_ix) = mi(step_dat, step_class(randperm(length(step_class))));
                    end
                    pval = (1+sum(miPerm >= thisMi)) / (analysisParams.statPerm + 1);
                    trajMiStatP{sess_ix, win_ix, fs_ix}(step_ix, 2:3) = [thisMi pval];
                end
            end %for each step
        end %for each feature set 
    end  %for each analysis window
    
    %% For baseline
    %Do for each featureset, using dimensionality reduction parameters.
    raster = getEpochRaster(cDat, baseAnaWin);
    baseFR = 1000*sum(raster > 0, 3)./sum(~isnan(raster), 3);
    thisMI = mi(baseFR, classId);
    for shuf_ix = 1:analysisParams.statPerm
        miPerm(shuf_ix) = mi(baseFR, classId(randperm(nTrials)));
    end
    pval = (1+sum(miPerm >= thisMi))/(analysisParams.statPerm+1);
    baseMiStatP(sess_ix, :) = [thisMi pval];
    
    fprintf(['Finished sessions ' num2str(sess_ix) ' of ' num2str(nSessions) '.\n']);
end

beep

save(fullfile(paths.results, 'mi_trajectories.mat'), 'analysisParams', 'trajMiStatP', 'baseMiStatP', 'featureNames');