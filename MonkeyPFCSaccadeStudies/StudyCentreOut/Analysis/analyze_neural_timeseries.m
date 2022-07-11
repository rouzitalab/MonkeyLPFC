% analyze_neural_timeseries
% Processes neural data one time-step through the trial at a time.

%% Paths and Constants
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths
sessions = my_sessions('CentreOut');
analysisParams = my_anaparams('timeseries');

%addpath(fullfile(paths.mlDevel, 'Helper'));
addpath(fullfile(paths.ml3rd, 'MIToolbox'));

%% Constants and other frequently accessed variables
stepSize = analysisParams.binWidth/2.5;  % 40 pcnt overlap
nAna = length(analysisParams.anaWins);
tsEdges = cell(1, nAna);
nSteps = nan(1, nAna);
for win_ix = 1:nAna
    ana_win = analysisParams.anaWins(win_ix);
    edgeStarts = ana_win.winEdges(1):stepSize:ana_win.winEdges(end)-1;
    edgeStops = edgeStarts + analysisParams.binWidth - 1;
    tsEdges{win_ix} = [edgeStarts' edgeStops'];
    while tsEdges{win_ix}(end, 2) > ana_win.winEdges(end) || tsEdges{win_ix}(end, 1) > ana_win.winEdges(end)
        tsEdges{win_ix} = tsEdges{win_ix}(1:end-1, :);
    end
    nSteps(win_ix) = size(tsEdges{win_ix}, 1);
end
clear win_ix ana_win edgeStarts edgeStops

nSessions = length(sessions);

%Preallocate cell arrays
timeseriesUnitMiP = cell(nSessions, nAna);  %Per-unit, per time-step MI and p-value
timeseriesUnitANP = cell(nSessions, nAna);  %Per-unit, per time-step F-stat and p-value, 
timeseriesUnitAN2 = cell(nSessions, nAna);  %Per-unit, all time-steps F-stat and p-value, using Mixed ANOVA
timeseriesUnitRMA = cell(nSessions, nAna);  %Per-unit, all time-steps F-stat and p-value, using RM ANOVA
timeseriesPopAN = cell(nSessions, nAna);    %All-units, per time-step F-stat and p-value, using manova1
timeseriesPopMi = cell(nSessions, nAna);    %All-units, per time-step MI and p-value
%Note, I could not figure out how to do all-units and all-time steps
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
    classId = [cDat.trial.newClass];
    [uqClasses, inClassIx] = unique(classId);
    nClasses = length(uqClasses);
    classStr = {cDat.trial(inClassIx).newClassStr};
    classDir = cat(1, cDat.trial(inClassIx).targPol); classDir = classDir(:, 1);
    clear inClassIx
    
    nUnits = size(cDat.trial(1).raster, 2);
    nTrials = length(cDat.trial);
    
    %% For each analysis window
    for win_ix = 1:nAna
        
        ana_win = analysisParams.anaWins(win_ix);
        cDat = trimTrials(cDat, ana_win);
        tVec = tsEdges{win_ix}(:,1) + diff(tsEdges{win_ix}, 1, 2)/2;
        
        %% Get the firing rate for each step for all units
        step_out_act = nan(nTrials, nUnits, nSteps(win_ix));
        for step_ix = 1:nSteps(win_ix)
            %Get the neural activity for this step for each trial
            tmp_anaWin = ana_win;
            tmp_anaWin.winEdges = tsEdges{win_ix}(step_ix, :);
            raster = getEpochRaster(cDat, tmp_anaWin);
            
            goodSampRast = squeeze(~any(isnan(raster), 2));
            nSampsPerTrial = sum(goodSampRast, 2);
            
            %If we want fRate
            nSpikesPerTrial = sum(raster, 3);
            FRate = 1000*bsxfun(@rdivide, nSpikesPerTrial, nSampsPerTrial);
            
            %keeping only trials with at least 70% the needed number of
            %samples
            sampCountBool = nSampsPerTrial > 0.7*size(raster, 3);
            
            %If we want spike counts, we need a common window
            firstSamp = find(~any(~goodSampRast(sampCountBool, :)), 1, 'first');
            lastSamp = find(~any(~goodSampRast(sampCountBool, :)), 1, 'last');
            nSpikesPerTrialCommon = sum(raster(:, :, firstSamp:lastSamp), 3);
            
            neur_act = sqrt(nSpikesPerTrialCommon);
            
            step_out_act(sampCountBool, :, step_ix) = neur_act(sampCountBool, :);
            
        end %step_ix
        clear step_ix tmp_anaWin raster goodSampRast nSampsPerTrial
        clear nSpikesPerTrial FRate firstSamp lastSamp nSpikesPerTrialCommon
        clear sampCountBool
        
        %% Do stats for each unit
        
        %Preallocate output
        timeseriesUnitANP{sess_ix, win_ix} = nan(nUnits, nSteps(win_ix), 2);
        timeseriesUnitMiP{sess_ix, win_ix} = nan(nUnits, nSteps(win_ix), 2);
        timeseriesUnitAN2{sess_ix, win_ix} = nan(nUnits, 2);
        timeseriesUnitRMA{sess_ix, win_ix} = nan(nUnits, 2);
        
        %%
        for unit_ix = 1:nUnits
            for step_ix = 1:nSteps(win_ix)
                %%
                this_step_act = step_out_act(:, unit_ix, step_ix);
                
                % Only keep trials that had data
                goodBool = ~isnan(this_step_act);
                this_step_act = this_step_act(goodBool);
                nGood = sum(goodBool);
            
                %Check that each class is well represented in this step.
                this_classIx = classId(goodBool)';
                [class_count, class_id] = hist(this_classIx, unique(this_classIx));
                
                if ~any(class_count < analysisParams.minTrialsPerGroup)
                    % Bin spike counts in (up to) 8 bins
                    tmp = this_step_act;
                    [uqAct, ~, bin_ix] = unique(tmp);
                    if length(uqAct) > 8
                        [~, binCenters] = hist(tmp(tmp>0), 7);
                        binEdges = [0 binCenters - diff([0 binCenters])./2 Inf];
                        [~, bin_ix] = histc(tmp, binEdges);
%                         bin_ix = ceil(8 * tiedrank(nSpikesPerTrialCommon) / length(nSpikesPerTrialCommon));
                    end
                    this_step_act = bin_ix;
                    
                    %Calculate ANOVA
                    [pval, tbl] = anova1(this_step_act, this_classIx, 'off');
                    timeseriesUnitANP{sess_ix, win_ix}(unit_ix, step_ix, :) = [tbl{2,5} pval];
                    
                    %Calculate mutual information.
                    thisMi = mi(this_step_act, this_classIx);
                    miPerm = nan(1, analysisParams.statPerm);
                    for shuf_ix = 1:analysisParams.statPerm
                        miPerm(shuf_ix) = mi(this_step_act, this_classIx(randperm(nGood)));
                    end
                    pval = (1+sum(miPerm >= thisMi)) / (analysisParams.statPerm + 1);
                    timeseriesUnitMiP{sess_ix, win_ix}(unit_ix, step_ix, :) = [thisMi pval];
                end
            end
            clear step_ix this_step_act goodBool nGood
            clear this_classIx class_count class_id
            clear pval tbl thisMi miPerm shuf_ix ix
            
            %% Do stats across all steps
            this_unit_act = squeeze(step_out_act(:, unit_ix, :));
            
            %Remove steps with any nan
            good_steps = ~any(isnan(this_unit_act));
            this_unit_act = this_unit_act(:, good_steps);
            this_nSteps = size(this_unit_act, 2);
            
            % ANOVA, within factor window, between factor class
            
            % Option A: 2-Way Mixed
            classFac = repmat(classId', 1, this_nSteps);  %Fixed effects
            winFac = repmat(tVec(good_steps)', nTrials, 1);  %Fixed effects
            trialFac = repmat((1:nTrials)', 1, this_nSteps);  %Random, nested in classFac
            myModel = [0 1 0; 0 0 1; 0 1 1];  %window, direction, window*direction
            data = this_unit_act(:);
            factors = [trialFac(:) winFac(:) classFac(:)];
            factors = factors(~isnan(data), :);
            [pvals,tbl,stats,anova_modl] = anovan(...
                data(~isnan(data)), {factors(:, 1), factors(:, 2), factors(:, 3)},...
                'model', myModel, 'sstype', 3,...
                'varnames', {'Trial' 'Window' 'Direction'},...
                'random', 1, 'nested', [0 0 1; 0 0 0; 0 0 0], 'display', 'off');
            timeseriesUnitAN2{sess_ix, win_ix}(unit_ix, :) = [tbl{4,6:7}];
            clear classFac winFac trialFac myModel data factors
            clear pvals tbl stats anova_modl
            
            % Option B: Two-way RM ANOVA
            varNames = cell(1, this_nSteps);
            for t_ix = 1:this_nSteps
                varNames{t_ix} = ['t' num2str(t_ix-1)];
            end
            varNames = cat(2, {'Direction'}, varNames);
            t = array2table([classId' this_unit_act],...
                'VariableNames', varNames);
            within_tbl = table((1:this_nSteps)','VariableNames',{'Window'});
            rm = fitrm(t,...
                [varNames{2} '-' varNames{end} ' ~ Direction'],...
                'WithinDesign', within_tbl);
            ranovatbl = ranova(rm);
            timeseriesUnitRMA{sess_ix, win_ix}(unit_ix, :) = [ranovatbl.F(2) ranovatbl.pValue(2)];
            clear varNames t_ix t within_tbl rm ranovatbl
            
            clear step_frate good_steps this_nSteps
        end; clear unit_ix
        
        %% For the population of neurons
        
        %Preallocate
        timeseriesPopAN{sess_ix, win_ix} = nan(nSteps(win_ix), 1+nClasses);
        timeseriesPopMi{sess_ix, win_ix} = nan(nSteps(win_ix), 2);
        
        %% Stats on each step
        for step_ix = 1:nSteps(win_ix)  %Last step is often problematic
            this_step_act = step_out_act(:, :, step_ix);
            
            % Remove any trials without this step
            good_trials = ~any(isnan(this_step_act), 2);
            this_step_act = this_step_act(good_trials, :);
            this_classId = classId(good_trials)';
            this_nTrials = sum(good_trials);
            
            % Remove any units with 0 spikes across all kept trials
            good_units = any(this_step_act ~= 0);
            this_step_act = this_step_act(:, good_units);
            
            %step_frate = uniqueCols(step_frate);
            
            if this_nTrials > 0 && size(this_step_act, 2) > 0 && ...
                    this_nTrials > size(this_step_act, 2)
                
                % manova1
                try
                    [d, pval] = manova1(this_step_act, this_classId);
                    timeseriesPopAN{sess_ix, win_ix}(step_ix, 1:length(pval)+1) = [d pval'];
                catch ME
                    if strcmpi(ME.identifier, 'stats:manova1:SingularSumSquares2')
                        fprintf('MANOVA1 on pop frate skipped for step %i\n', step_ix);
                    else
                        rethrow(ME)
                    end
                end

                % mi - requires lots of memory
                thisMi = mi(this_classId, this_step_act);
                miPerm = nan(analysisParams.statPerm, 1);
                for shuf_ix = 1:analysisParams.statPerm
                    miPerm(shuf_ix) = mi(this_classId(randperm(this_nTrials)), this_step_act);
                end
                pval = (1+sum(miPerm >= thisMi))/(analysisParams.statPerm+1);
                timeseriesPopMi{sess_ix, win_ix}(step_ix, :) = [thisMi pval];
            
            end
        end
        clear step_ix step_frate d pval
        clear thisMi miPerm shuf_ix
        
        %% TODO: Stats on entire timeseries
        % How?
        
    end
    clear unit_ix
    fprintf(['Finished sessions ' num2str(sess_ix) ' of ' num2str(nSessions) '.\n']);
end

beep

save(fullfile(paths.results, 'mi_timeseries.mat'),...
    'analysisParams',...
    'timeseriesUnitMiP', 'timeseriesUnitANP', ...
    'timeseriesUnitAN2', 'timeseriesUnitRMA', ...
    'timeseriesPopAN', 'timeseriesPopMi');