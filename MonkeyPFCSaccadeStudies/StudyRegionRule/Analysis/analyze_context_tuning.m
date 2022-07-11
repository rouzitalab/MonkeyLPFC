% analyze_context_tuning.m
% Looks at each unit in each session to see if its activity is different
% between rules (contexts) for behaviours to the same target. This applies
% only to saccades to targets in the corners.
% For each corner:
% Firing rate is calculated in the visual period (targets on-screen, but
% not the cue) and the cue period (rule known).
% A paired t-test is done for each period (rule A vs rule B)
% An ANOVA is done to see if firing rate differences in cue period are
% greater than firing rate differences in visual period.
% All firing rates are normalized by within-block baseline.

%% Paths and Constants
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths
sessions = my_sessions('RegionRule');
analysisParams = my_anaparams('ContextTuning');

%addpath(fullfile(paths.mlcbb, 'Helper'));
%addpath(fullfile(paths.ml3rd, 'MIToolbox'));

%% Constants
nAna = length(analysisParams.anaWins);  % Number of analysis windows.
nCorners = 4;  % Targets in the corners
nRules = 2;  % 2 rules per corner

%% Pre-allocate variables we will save.
nSessions = length(sessions);
baselineFR = cell(nSessions, 1);
tuningStatP = cell(nSessions, 1);
tuningCurve = cell(nSessions, 1);
timeseriesStatP = cell(nSessions, 1);
trialBoolOut = cell(nSessions, 1);
fullInvalidOut = cell(nSessions, 1);
outFRate = cell(nSessions, 1);

% quick init for debugging
sess_ix = 1;
unit_ix = 17;

%%
for sess_ix = 1:nSessions
    %%
    this_sess = sessions(sess_ix);
    
    %% Load and process ptbmat
    ptb = load(fullfile(paths.preprocessed, 'ptb', this_sess.ptbfname), '-mat');
    [ptb.eData.trial, ptb.params] = getTrialStimInfo(ptb.eData.trial, ptb.params);
    ptb.eData.trial = getTrialBehavResult(ptb.eData.trial, ptb.params, analysisParams);
    [IsInCondition] = getCornerConditionMembership(ptb.eData.trial);
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
    IsInCondition = IsInCondition(triageBool, :, :, :);
    trialBoolOut{sess_ix} = triageBool;
    cDat.trial = cDat.trial(triageBool); clear trialBool
    nTrials = length(cDat.trial);
    
    % Eliminate units that fire < 1Hz
    % Also verifies this to be true per-block
    % Adds cDat.invalidUnits.
    % If thisParams.useUnits == 'merged' then it merges all units per
    % channel into one unit per channel.
    %53 units including noise.
    %19 sorted units after triage.
    [cDat, fullInvalidOut{sess_ix}] = triageUnits(cDat, analysisParams);
    nUnits = size(cDat.trial(1).raster, 2);
    
    
    %% Pre-allocate output variables
    outFRate{sess_ix} = nan(nUnits, nAna, nCorners, nRules);
    
    unit_ix = 1; win_ix = 1; corner_ix = 1; rule_ix = 1; % Quick init for debugging.
    %% For each unit in this session
    for unit_ix = 1:nUnits
        %% For each analysis window,
        for win_ix = 1:nAna
            %%
            ana_win = analysisParams.anaWins(win_ix);
            raster = getRasterForWindow(cDat, ana_win);
            
            % Save the baseline firing rate
            if win_ix == 1
                baselineFR{sess_ix}(unit_ix) = 1000*sum(nansum(raster, 2))/sum(sum(~isnan(raster)));
            end
            
            % Some vectors we'll use to slice up the raster.
            plotXVec = ana_win.plotX(1):ana_win.plotX(2);
            xEdges = ana_win.plotX(1)-1:analysisParams.binWidth:ana_win.plotX(2);
            
            %% Per-corner & Rule
            for corner_ix = 1:nCorners
                for rule_ix = 1:nRules
                    %%
                    % Trials with this class
                    trBool = strcmpi({cDat.trial.newClassStr}, classStr{cl_ix});
                    
                    %Prepare Raster: scatter(spkTimes,spkTrial)
                    tmpRaster = raster(trBool, :);
                    [spkTrial, ti] = find(tmpRaster > 0);  % trial index and time index of each spike.
                    spkTimes = plotXVec(ti);
                    
                    %Prepare PSTH: bar(xEdges, rtSpkBin)
                    nSpkBin = histc(spkTimes, xEdges);  % Bin the spike times according to xEdges
                    [~, tsamp] = find(~isnan(tmpRaster));  % indices of ~nan samples
                    sampTimes = plotXVec(tsamp);  % times of ~nan samples
                    nSampBin = histc(sampTimes, xEdges);  % number of samples per bin
                    rtSpkBin = 1000*nSpkBin./nSampBin;  % spike rate (spikes/sec)
                    dirFRate = 1000*sum(sum(tmpRaster > 0))/sum(sum(~isnan(tmpRaster)));%Total spiking rate for baseline normalization and tuning curves.
                    rtSpkBin = rtSpkBin * (1/tuningCurve{sess_ix}(unit_ix, 1, cl_ix)) * (maxNTrials/3) * (1/5);  % Scale by baseline. Plot max = 5x baseline.
                    
                    %Save per-direction fRate into tuningCurve
                    tuningCurve{sess_ix}(unit_ix, win_ix, cl_ix) = dirFRate;
                    
                    
                end
            end
            %%
            clear plotXVec xEdges cl_ix trBool tmpRaster spkTrial ti spkTimes
            clear nSpkBin tsamp sampTimes nSampBin rtSpkBin dirFRate
            
            
            
            %% Do stats on (unbinned) spike rate.
            trFRate = 1000*sum(raster > 0, 2)./sum(~isnan(raster), 2);
            nTrials = length(trFRate);
            
            % MI with permutations to get p
            thisMi = mi(trFRate,classId');
            miPerm = nan(analysisParams.statPerm, 1);
            for shuf_ix = 1:analysisParams.statPerm
                miPerm(shuf_ix) = mi(trFRate, classId(randperm(nTrials))');
            end
            pval = (1+sum(miPerm >= thisMi))/(analysisParams.statPerm+1);
            tuningStatP{sess_ix}(unit_ix, win_ix, :) = [thisMi pval];
            
            clear trFRate nTrials thisMi miPerm shuf_ix ix pval
            
            % Alternatives
            %
            % kruskalwallis (non-parametric one-way ANOVA)
            % [pvals, tbl, stats] = kruskalwallis(trFRate, classIx, 'off');
            % tuningStatP{sess_ix}(unit_ix, win_ix-1, :) = [tbl{2,5} pvals];
            %
            % 2-way mixed ANOVA, within factor window, between factor class
            % Option A:
            % classFac = [classIx';classIx'];  % fixed
            % winFac = [ones(nTrials,1);2*ones(nTrials,1)];  % fixed
            % trialFac = [(1:nTrials)';(1:nTrials)'];  % random, nested in classFac
            % myModel = [0 1 0; 0 0 1; 0 1 1];
            % [pvals,tbl,~,~] = anovan([baselineFRate;trFRate], {trialFac winFac classFac},...
            %     'model', myModel, 'sstype', 3, 'varnames', {'trialInd' 'window' 'direction'},...
            %     'random', 1, 'nested', [0 0 1; 0 0 0; 0 0 0], 'display', 'off');
            % Option B:
            % t = table(classIx', baselineFRate, trFRate,...
            %     'VariableNames',{'Direction','t0','t1'});
            % rm = fitrm(t,'t0-t1 ~ Direction','WithinDesign', [1 2]');
            % ranovatbl = ranova(rm);
            % tuningStatP{sess_ix}(unit_ix, win_ix-1, :) = [tbl{4, 6} pvals(3)];
        end
        clear win_ix
        
        
        clear sph totSpkRate
        
    end
    clear unit_ix
    
    %     %% Calculate the proportion of units with significant tuning
    %     for win_ix = 1:size(tuningKWChiSqP, 2)
    %         uBool = tuningKWChiSqP(:, win_ix, 2) < 0.01;
    %         pTuned(sess_ix, win_ix, :) = hist(tuningPrefDir(uBool, win_ix), 1:8) ./ size(tuningPrefDir, 1);
    %     end
    %     nUnits(sess_ix) = size(tuningPrefDir, 1);
    %     clear win_ix uBool
    
    %     %% Save data for DataHigh
    %     %D(trial_ix).
    %     %   data (nUnits x n_1msec_bins raster)
    %     % [the rest are optional]
    %     %   type ('traj' or 'state')
    %     %   epochStarts (1xnEpochs; indices of the starts of each epoch)
    %     %   epochColors (nEpochsxRGB)
    %     %   condition ('conditionName')
    %
    %     addpath(fullfile(paths.ml3rd, 'DataHigh1.1'));
    %     thisParams.timeLockEvent = 'target onset';  %'target onset' 'fixation offset' 'saccade onset'
    %     thisParams.winEdges = [-100 2000];
    %     thisParams.avoidEvent = 'saccade onset';
    %     thisParams.avoidWin = [-50 Inf];
    %
    %     cDat = trimTrials(cDat, thisParams);% Identify eventTime, avoidTime, a tVec, and timeBool for each trial.
    %     D = cDat2DataHigh(cDat.trial);
    %     DataHigh(D, 'DimReduce');
    %     rmpath(fullfile(paths.ml3rd, 'DataHigh1.1'));
    fprintf(['Finished ' num2str(sess_ix) ' of ' num2str(nSessions) '.\n']);
end
beep

%% Reshape data
sessId = [];
for sess_ix = 1:nSessions
    sessId = [sessId; sess_ix*ones(length(baselineFR{sess_ix}), 1)];
end
clear sess_ix

baselineFR = cat(1, baselineFR{:});  % In case I want to further trim units.
tuningCurve = cat(1, tuningCurve{:});
tuningStatP = cat(1, tuningStatP{:});
timeseriesStatP = cat(1, timeseriesStatP{:});

%%
save(fullfile(paths.results, 'tuning.mat'), 'trialBoolOut', 'fullInvalidOut', 'baselineFR', 'tuningCurve', 'tuningStatP', 'sessId', 'analysisParams', 'sessions');