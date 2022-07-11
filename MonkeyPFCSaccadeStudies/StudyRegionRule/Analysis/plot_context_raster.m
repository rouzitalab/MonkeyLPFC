% neural_tuning
% A script to load preprocessed behavioural and neural data, consolidate,
% then describe each neuron by its directional selectivity.
%TODO: Move the exemplar plotting to a separate script.

%% Paths and Constants
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths
sessions = my_sessions('RegionRule');
analysisParams = my_anaparams('tuning');

nAna = length(analysisParams.anaWins);  % Number of analysis windows.
% analysisParams.statPerm = 1;  % Override for quick plotting.

%addpath(fullfile(paths.mlcbb, 'Helper'));
addpath(fullfile(paths.ml3rd, 'MIToolbox'));

%% Specify plotting details.
%Map the corner to plot location
xoff = 2/100; yoff = 2/100; w = 0.45; h = 0.45;
plotP.pos = [...
    0.5+xoff, 0.5+yoff, w, h;...  % UR
    0.5+xoff,     yoff, w, h;...  % DR
        xoff,     yoff, w, h;...  % DL
        xoff, 0.5+yoff, w, h];    % UL
clear xoff yoff w h tmp_ix
plotP.lineStyle = {'-', '--'};
plotP.psth_scale = 1/3;
plotP.psth_size = 1/3;
plotP.exemplars = [2 15; 5 25];

%% Pre-allocate variables we will save.
nSessions = length(sessions);
baselineFR = cell(nSessions, 1);
tuningStatP = cell(nSessions, 1);
tuningCurve = cell(nSessions, 1);
timeseriesStatP = cell(nSessions, 1);
trialBoolOut = cell(nSessions, 1);
fullInvalidOut = cell(nSessions, 1);
    
% quick init for debugging
sess_ix = 4;
unit_ix = 25;

%%
for sess_ix = 1:nSessions
    %%
    this_sess = sessions(sess_ix);
    
    %% Load and process ptbmat
    ptb = load(fullfile(paths.preprocessed, 'ptb', this_sess.ptbfname), '-mat');
    [ptb.eData.trial, ptb.params] = getTrialStimInfo(ptb.eData.trial, ptb.params);
    ptb.eData.trial = getTrialBehavResult(ptb.eData.trial, ptb.params, analysisParams);
    [IsInCondition] = getCornerConditionMembership(ptb.eData.trial);
    ptb.eData.trial = getNewClass(ptb.eData.trial, struct('classifyTarg', 'targClass'));
    
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
    IsInCondition = IsInCondition(trialBool, :, :, :);
    trialBoolOut{sess_ix} = trialBool;
    cDat.trial = cDat.trial(trialBool); clear trialBool
    
    % Eliminate units that fire < 1Hz or > 100 Hz.
    % Adds cDat.invalidUnits.
    % If thisParams.useUnits == 'merged' then it merges all units per
    % channel into one unit per channel.
    %53 units including noise.
    %19 sorted units after triage.
    [cDat, fullInvalidOut{sess_ix}] = triageUnits(cDat, analysisParams);
    
    
    %% Pre-allocate output variables
    nUnits = size(cDat.trial(1).raster, 2);
    nCorners = 4;
    nRules = 2;
    nTrials = sum(sum(any(IsInCondition, 4)), 3);
    
    %% For each unit in this session
    for unit_ix = 1:nUnits
        %%
        win_ix = 1;

        %% Prepare figure
        figure('Name', ['Session ' num2str(sess_ix) ', Unit ' num2str(unit_ix)]);
        
        sph = nan(1, 4);  %+1 for tuning curve in centre.
        for sp_ix = 1:4
            sph(sp_ix) = subplot('Position', plotP.pos(sp_ix, :));
            box off
            set(gca,'Color', 'none')
            set(gca, 'XColor', 'none')
            set(gca, 'YColor', 'none')
            hold on
        end
        clear sp_ix
        
        tuningCurve = nan(nAna, nCorners, nRules);
        %% For each analysis window,
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
            
            % Some vectors we'll use to slice up the raster.
            plotXVec = ana_win.plotX(1):ana_win.plotX(2);
            xEdges = ana_win.plotX(1)-1:analysisParams.binWidth:ana_win.plotX(2);
            
            %% Per-corner&rule: raster, PSTH, overall fRate for tuningCurve.
            rtSpkBinOut = nan(length(xEdges), 4, 2);
            for corner_ix = 1:4
                rule_trials = [0;squeeze(sum(any(IsInCondition(:, corner_ix, :, :), 4)))];
                for rule_ix = 1:2
                    trBool = squeeze(any(IsInCondition(:, corner_ix, rule_ix, :), 4));
                    
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
                    tuningCurve(win_ix, corner_ix, rule_ix) = dirFRate;
                    rtSpkBinOut(:, corner_ix, rule_ix) = rtSpkBin * (1/tuningCurve(1, corner_ix, rule_ix)) * (max(nTrials)*plotP.psth_size) * plotP.psth_scale;  % Scale by baseline. Plot max = 5x baseline.
                    
                    plot(sph(corner_ix), spkTimes, ones(1, length(spkTimes))*rule_trials(rule_ix), 'k-')
                    scatter(sph(corner_ix), spkTimes, spkTrial+rule_trials(rule_ix), [ana_win.colour '.']);
                end
            end
            %
            ylims = [min(rtSpkBinOut(:)) max(rtSpkBinOut(:))];
            for corner_ix = 1:4
                for rule_ix = 1:2
                    if rule_ix == 1
                        fc = [0.8 0.8 0.8];
                    else
                        fc = 'none';
                    end
                    h = bar(sph(corner_ix), xEdges(1:end-1) + analysisParams.binWidth/2,...
                        -1*squeeze(rtSpkBinOut(1:end-1, corner_ix, rule_ix)),...
                        1, 'FaceColor', fc, 'EdgeColor', ana_win.colour, 'LineWidth', 1);
                    
                end
                plot(sph(corner_ix), [ana_win.vBar ana_win.vBar], [-max(nTrials)*plotP.psth_size max(nTrials)], 'k--');
                set(sph(corner_ix), 'YDir', 'reverse')
                set(sph(corner_ix), 'YLim', [-max(nTrials)*plotP.psth_size max(nTrials)]);
            end
            %%
            clear plotXVec xEdges cl_ix trBool tmpRaster spkTrial ti spkTimes
            clear nSpkBin tsamp sampTimes nSampBin rtSpkBin dirFRate
            
        end
        clear win_ix
        
        
        clear sph totSpkRate
        
    end
    clear unit_ix
    fprintf(['Finished ' num2str(sess_ix) ' of ' num2str(nSessions) '.\n']);
end