% neural_tuning
% Describe each neuron by its directional selectivity.

%% Paths and Constants
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths
sessions = my_sessions('CentreOut');
analysisParams = my_anaparams('tuning');
fullTrialAnaParams = my_anaparams('classification');
fullTrialAnaParams.anaWins = fullTrialAnaParams.anaWins(...
    strcmpi({fullTrialAnaParams.anaWins.name}, 'fullTrial'));

MIN_REL_WINDOW_SIZE = 0.7;  % Amount of window that must have good samples.
STAT_VAR = 'bin_ix';  % What variable should we do stats on?
                      % Options: bin_ix, FRate, nSpikesPerTrialCommon
                      % sqrt(FRate), sqrt(nSpikesPerTrialCommon)

%addpath(fullfile(paths.mlcbb, 'Helper'));
addpath(fullfile(paths.ml3rd, 'MIToolbox'));
addpath(fullfile(paths.ml3rd, 'RossMI'));

%% Constants
nAna = length(analysisParams.anaWins);  % Number of analysis windows.
% analysisParams.statPerm = 1;  % Override for quick plotting.

%% Pre-allocate variables we will save.
nSessions = length(sessions);
tuningFRate = cell(nSessions, 1);       % Helps know why a unit was discarded
tuningCountVar = cell(nSessions, 1);    % Keep track of mean spike count and variance
tuningStatAN = cell(nSessions, 1);      % ANOVA F-stat and p-value
tuningStatANMC = cell(nSessions, 1);    % t-test for 1 vs rest
tuningStatMI = cell(nSessions, 1);      % Mutual Inf. and perm. p-value
tuningStatKW = cell(nSessions, 1);      % Kruskal-Wallis stat and p-value
tuningStatPopAN = cell(nSessions, 1);   % Pop ANOVA F-stat and p-value
tuningStatPopMI = cell(nSessions, 1);   % Pop MI and p-value
tuningCurve = cell(nSessions, 1);       % Firing rates
trialBoolOut = cell(nSessions, 1);      % Which trials were kept
fullInvalidOut = cell(nSessions, 1);    % Which units were kept
sortQuality = cell(nSessions, 1);       % The quality of the sort from mksort
meanFRate = cell(nSessions, 1);         % The total FRate across all samples and trials

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
    
    %plotFigure1(ptb.eData.trial(174), cDat.trial(174))
    
    % Eliminate trials
    % - without saccades
    % - without timeLockEvent
    % - without appropriate newOutcomeCodes
    % - without minTrialsPerGroup
    trialBool = triageTrials(cDat, analysisParams);
    trialBoolOut{sess_ix} = trialBool;
    
    cDat.trial = cDat.trial(trialBool); clear trialBool
    
    % Select units
    % Adds cDat.invalidUnits. See (global) DEF.unitState
    % If analysisParams.useUnits == 'merged' then it merges all units per
    % channel into one unit per channel (multiunit). (['all'], 'sorted')
    % If analysisParams.minAllowedBaselineFR, it eliminates neurons that
    % don't meet this minimum firing rate threshold.
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
    % Epoch firing rates. If it == 0, then no spikes in any trial.
    tuningFRate{sess_ix} = nan(nUnits, nAna);
    % per ana_win tuning curve; from this we can get preferred direction.
    tuningCurve{sess_ix} = nan(nUnits, nAna, nClasses);
    tuningCountVar{sess_ix} = nan(nUnits, nAna, 2);
    
    tuningStatAN{sess_ix} = nan(nUnits, nAna, 2);% F and p
    tuningStatANMC{sess_ix} = nan(nUnits, nAna, nClasses);
    tuningStatMI{sess_ix} = nan(nUnits, nAna, 2);% chisq & p
    tuningStatKW{sess_ix} = nan(nUnits, nAna, 2);
    tuningStatPopAN{sess_ix} = nan(nAna, nClasses);
    tuningStatPopMI{sess_ix} = nan(nAna, 2);
    
    % Some other epoch-independent variables
    % quality of each unit's sort (1-4)
    sortQuality{sess_ix} = cDat.sortQuality;
    % total firing rate across all trials (epoch independent)
    raster = getEpochRaster(cDat, fullTrialAnaParams.anaWins);
    nsamps = sum(sum(~isnan(raster), 3));
    spkCounts = nansum(nansum(raster, 3));
    meanFRate{sess_ix} = 1000*spkCounts./nsamps;
    
    %% For each unit in this session
    fprintf('Analyzing tuning in %i units...\n', nUnits);
    for unit_ix = 1:nUnits
        %% For each analysis window,
        % Do raster plots (maybe), get firing rate, chisq, p, tuning curve
        for win_ix = 1:nAna
            %%
            ana_win = analysisParams.anaWins(win_ix);
            raster = getEpochRaster(cDat, ana_win);  % trials x units x samples
            raster = squeeze(raster(:, unit_ix, :));  % trials x samples
            
            % Extract some simple variables
            goodSampRast = ~isnan(raster);
            nSampsPerTrial = sum(goodSampRast, 2);
            nSpikesPerTrial = sum(raster > 0, 2);
            
            tuningCountVar{sess_ix}(unit_ix, win_ix, :) =...
                [mean(nSpikesPerTrial) var(nSpikesPerTrial)];

            % Save the firing rate for all trials in epoch
            % Epochs with firingRate == 0 will be excluded from some
            % analyses.
            tuningFRate{sess_ix}(unit_ix, win_ix) =...
                1000*sum(nSpikesPerTrial)/sum(nSampsPerTrial);
            
            %% Per-direction: raster, PSTH, overall fRate for tuningCurve.
            
            % trials containing good samples for at least 70% the window
            sampCountBool = nSampsPerTrial >...
                floor(MIN_REL_WINDOW_SIZE*size(raster, 2));
            % for each direction
            for cl_ix = 1:nClasses
                % Keep trials with this class and a sufficiently long win
                trBool = sampCountBool &...
                    strcmpi({cDat.trial.newClassStr}, classStr{cl_ix})';
                
                %Save per-direction fRate into tuningCurve
                % Some trials have incomplete epochs, so only consider good
                % samples.
                if any(trBool)
                    tuningCurve{sess_ix}(unit_ix, win_ix, cl_ix) = ...
                        1000*sum(nSpikesPerTrial(trBool))/sum(nSampsPerTrial(trBool));
                end
            end
            clear cl_ix trBool

            %% Do stats on spiking activity
            % Only if there is non-zero activity in this epoch.
            if tuningFRate{sess_ix}(unit_ix, win_ix) > 0
                
                % Identify class ids for all trials in epoch
                this_classId = classId(sampCountBool)';  % Target locations
                nTrials = length(this_classId);
                
                % Calculate some variables we might use for stats
                
                % Spike counts - variable length of epochs across trials
                % require normalization: FRate
                nSpikesPerTrial = nSpikesPerTrial(sampCountBool);
                nSampsPerTrial = nSampsPerTrial(sampCountBool);
                FRate = 1000*nSpikesPerTrial./nSampsPerTrial;  %FRate
                
                % Spike counts - common length trials
                firstSamp = find(~any(~goodSampRast(sampCountBool, :)), 1, 'first');
                lastSamp = find(~any(~goodSampRast(sampCountBool, :)), 1, 'last');
                nSpikesPerTrialCommon = sum(raster(sampCountBool, firstSamp:lastSamp), 2);
                
                % Binned spike counts (or binned FRate). Up to 8 bins.
                tmp = sqrt(nSpikesPerTrialCommon);  %sqrt(FRate)
                [uqAct, ~, bin_ix] = unique(tmp);
                if length(uqAct) > 8
                    [~, binCenters] = hist(tmp(tmp>0), 7);
                    binEdges = [0 binCenters - diff([0 binCenters])./2 Inf];
                    [~, bin_ix] = histc(tmp, binEdges);
                    %                 bin_ix = ceil(8 * tiedrank(nSpikesPerTrialCommon) / length(nSpikesPerTrialCommon));
                end
                
                neurAct = eval(STAT_VAR);
                
                %ANOVA on firing rate bin with factor class
                [pval, tbl, stats] = anova1(neurAct, this_classId, 'off');
                tuningStatAN{sess_ix}(unit_ix, win_ix, :) = [tbl{2,5} pval];
                %[c,m,h] = multcompare(stats, 'CType', 'bonferroni');
                
                %Compare each class vs all other classes
                for cl_ix = 1:nClasses
                    this_bool = this_classId==cl_ix;
                    [h, p, ci, stats] = ttest2(neurAct(this_bool), neurAct(~this_bool));
                    tuningStatANMC{sess_ix}(unit_ix, win_ix, cl_ix) = ...
                        p*sign(mean(neurAct(this_bool))-mean(neurAct(~this_bool)));
                end
                
                % MI with permutations to get p
                %             [thisMi, ~] = discrete_continuous_info(this_classId, neurAct, 3, 2);
                %             miPerm = nan(analysisParams.statPerm, 1);
                %             for shuf_ix = 1:analysisParams.statPerm
                %                 [miPerm(shuf_ix), ~] = discrete_continuous_info(...
                %                     this_classId(randperm(nTrials)), neurAct, 3, 2);
                %             end
                
                thisMi = mi(this_classId, neurAct);
                miPerm = nan(analysisParams.statPerm, 1);
                for shuf_ix = 1:analysisParams.statPerm
                    miPerm(shuf_ix) = mi(this_classId(randperm(nTrials)), neurAct);
                end
                
                pval = (1+sum(miPerm >= thisMi))/(analysisParams.statPerm+1);
                tuningStatMI{sess_ix}(unit_ix, win_ix, :) = [thisMi pval];
                
                % kruskalwallis (non-parametric one-way ANOVA)
                [pval, tbl] = kruskalwallis(neurAct, this_classId, 'off');
                tuningStatKW{sess_ix}(unit_ix, win_ix, :) = [tbl{2,5} pval];
                
            else
                tuningStatAN{sess_ix}(unit_ix, win_ix, :) = [nan nan];
                tuningStatANMC{sess_ix}(unit_ix, win_ix, :) = nan;
                tuningStatMI{sess_ix}(unit_ix, win_ix, :) = [nan nan];
                tuningStatKW{sess_ix}(unit_ix, win_ix, :) = [nan nan];
            end
            
            clear goodSampRast nSampsPerTrial nSpikesPerTrial sampCountBool
            clear trFRate nTrials thisMi miPerm shuf_ix pval tbl stats this_classId
        end
        clear win_ix
        clear sph totSpkRate
        fprintf('.');
    end
    clear unit_ix
    fprintf('Done.\n');

    %% Population analysis
    
    % Only keep units that have non-zero activity in each epoch.
    nonzerobool = ~any(tuningFRate{sess_ix} == 0, 2);
    
    for win_ix = 1:nAna
        
        ana_win = analysisParams.anaWins(win_ix);
        raster = getEpochRaster(cDat, ana_win);
        raster = raster(:, nonzerobool, :);
        
        % Only keep trials that have at least MIN_REL_WINDOW_SIZE of the epoch
        nGoodSamps = sum(~isnan(raster), 3);
        trBool = nGoodSamps(:,1) > floor(MIN_REL_WINDOW_SIZE*size(raster, 3));
        raster = raster(trBool, :, :);
        nGoodSamps = nGoodSamps(trBool, :);
        this_classId = classId(trBool)';
        %nTrials = size(this_classId, 1);
        
        % Firing rate (using variable-length trials)
        FRate = 1000*sum(raster > 0, 3)./nGoodSamps;
        
        % Spike counts (trim trials to common-length)
        goodSamps = find(~squeeze(any(isnan(raster(:, 1, :)), 1)));
        nSpikesPerTrialCommon = sum(raster(:, :, goodSamps(1):goodSamps(end)) > 0, 3);
        
        % Bin into up to 8 bins
        tmp = sqrt(nSpikesPerTrialCommon);  % sqrt(FRate)
        [uqAct, ~, bin_ix] = unique(tmp);
        if length(uqAct) > 8
            [~, binCenters] = hist(tmp(tmp>0), 7);
            binEdges = [0 binCenters - diff([0 binCenters])./2 Inf];
            [~, bin_ix] = histc(tmp, binEdges);
%             bin_ix = ceil(8 * tiedrank(nSpikesPerTrialCommon) / length(nSpikesPerTrialCommon));
        end
        
        % Get preferred output
        neurAct = eval(STAT_VAR);
        
        % Remove any columns (neurons) with no variability.
        nonzerovarbool = var(neurAct) ~= 0;
        
        %Population analysis using manova1
        [d, pval] = manova1(neurAct(:, nonzerovarbool), this_classId);
        tuningStatPopAN{sess_ix}(win_ix, :) = [d pval'];
        
        %Population analysis using multivariate MI
        
%         [thisMi, ~] = discrete_continuous_info(this_classId, neurAct);
%         miPerm = nan(analysisParams.statPerm, 1);
%         for shuf_ix = 1:analysisParams.statPerm
%             [miPerm(shuf_ix), ~] = discrete_continuous_info_fast(...
%                 this_classId(randperm(nTrials)), neurAct);
%         end
            
%         thisMi = mi(this_classId, nSpikesPerTrialCommon);
%         miPerm = nan(analysisParams.statPerm, 1);
%         for shuf_ix = 1:analysisParams.statPerm
%             miPerm(shuf_ix) = mi(this_classId(randperm(nTrials)), nSpikesPerTrialCommon);
%         end
%         
%         pval = (1+sum(miPerm >= thisMi))/(analysisParams.statPerm+1);
%         tuningStatPopMI{sess_ix}(win_ix, :) = [thisMi pval];
    end
    
    clear win_ix ana_win tVec
    clear trFRate nTrials d pval stats thisMi miPerm shuf_ix
    
    fprintf(['Finished session ' num2str(sess_ix) ' of ' num2str(nSessions) '.\n']);
end
beep

%%
save(fullfile(paths.results, 'tuning_epochs.mat'),...
    'trialBoolOut', 'fullInvalidOut', 'tuningFRate', 'tuningCountVar',...
    'tuningCurve', 'sortQuality', 'meanFRate',...
    'tuningStatAN', 'tuningStatANMC', 'tuningStatMI', 'tuningStatKW',...
    'tuningStatPopAN', 'tuningStatPopMI',...
    'analysisParams', 'sessions');