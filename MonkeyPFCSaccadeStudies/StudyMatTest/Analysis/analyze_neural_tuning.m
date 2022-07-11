% neural_tuning
% A script to load preprocessed behavioural and neural data, consolidate,
% then describe each neuron by its directional selectivity.
%TODO: Move the exemplar plotting to a separate script.

%% Paths and Constants
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths
sessions = my_sessions('Contains_Inverse');
analysisParams = my_anaparams('tuning');

%addpath(fullfile(paths.mlcbb, 'Helper'));
addpath(fullfile(paths.ml3rd, 'MIToolbox'));

%% Constants
nAna = length(analysisParams.anaWins);  % Number of analysis windows.
% analysisParams.statPerm = 1;  % Override for quick plotting.

%% Pre-allocate variables we will save.
nSessions = length(sessions);
baselineFR = cell(nSessions, 1);
tuningStatP = cell(nSessions, 1);
tuningCurve = cell(nSessions, 1);
timeseriesStatP = cell(nSessions, 1);
trialBoolOut = cell(nSessions, 1);
fullInvalidOut = cell(nSessions, 1);
    
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
    trialBoolOut{sess_ix} = trialBool;
    cDat.trial = cDat.trial(trialBool); clear trialBool
    
    % Eliminate units that fire < 1Hz or > 100 Hz.
    % Adds cDat.invalidUnits.
    % If thisParams.useUnits == 'merged' then it merges all units per
    % channel into one unit per channel.
    %53 units including noise.
    %19 sorted units after triage.
    [cDat, fullInvalidOut{sess_ix}] = triageUnits(cDat, analysisParams);
    
    %% Easy-access, map classes -> directions
    classId = [cDat.trial.newClass];
    [uqClasses, inClassIx] = unique(classId);
    nClasses = length(uqClasses);
    classStr = {cDat.trial(inClassIx).newClassStr};
    classDir = cat(1, cDat.trial(inClassIx).targPol); classDir = classDir(:, 1);
    clear inClassIx
    
    nUnits = size(cDat.trial(1).raster, 2);
    maxNTrials = max(hist(classId, uqClasses));  % Max number of trials per direction to set y-axis.
    
    %% Pre-allocate output variables
    % Baseline firing rates (to know why a unit was discarded)
    baselineFR{sess_ix} = nan(nUnits, 1);
    % per ana_win tuning curve; from this we can get preferred direction.
    tuningCurve{sess_ix} = nan(nUnits, nAna, nClasses);
    % per ana_win chisq & p
    tuningStatP{sess_ix} = nan(nUnits, nAna, 2);
    
    %% For each unit in this session
    for unit_ix = 1:nUnits
        %%
        
        
        
        win_ix = 1;
        %% For each analysis window,
        % Do raster plots (maybe), get firing rate, chisq, p, tuning curve
        for win_ix = 1:nAna
            %%
            ana_win = analysisParams.anaWins(win_ix);
            
            % Identify eventTime, avoidTime, a tVec, and timeBool for each trial.
            cDat = trimTrials(cDat, ana_win);
            
            %Plug each trial's raster into a common raster.
            tVec = ana_win.winEdges(1):ana_win.winEdges(2);
            raster = nan(length(cDat.trial), length(tVec));
            for tr_ix = 1:length(cDat.trial)
                rBool = cDat.trial(tr_ix).timeBool;  % this trial's in-window samples
                trTVec = cDat.trial(tr_ix).tVec(rBool);  % the times of this trial's in-window samples
                lBool = tVec >= trTVec(1) & tVec <= trTVec(end);  % the samples in the common raster
                raster(tr_ix, lBool) = cDat.trial(tr_ix).raster(rBool, unit_ix);
            end
            clear tVec tr_ix rBool trTVec lBool
            
            % Save the baseline firing rate
            if win_ix == 1
                baselineFR{sess_ix}(unit_ix) = 1000*sum(nansum(raster, 2))/sum(sum(~isnan(raster)));
            end
            
            % Some vectors we'll use to slice up the raster.
            plotXVec = ana_win.plotX(1):ana_win.plotX(2);
            xEdges = ana_win.plotX(1)-1:analysisParams.binWidth:ana_win.plotX(2);
            
            cl_ix = 1;
            %% Per-direction: raster, PSTH, overall fRate for tuningCurve.
            for cl_ix = 1:nClasses
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