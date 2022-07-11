% This script looks at the neural trajectories obtained using GPFA and
% manova1 on different epochs in a little more detail.

%% Paths
% addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
% my_consts;
% my_paths;
% global paths DEF
% 
% addpath(genpath(fullfile(paths.ml3rd, 'DataHigh')));
% addpath(fullfile(paths.ml3rd, 'dPCA', 'matlab'));

%% Constants

sessions = my_sessions('CentreOut');
nSessions = length(sessions);

analysisParams = my_anaparams('trajOnly');
trajAnaWin = analysisParams.anaWins;

epochAnaParams = my_anaparams('tuning');
epochAnaWin.base = epochAnaParams.anaWins(strcmpi({epochAnaParams.anaWins.name}, 'baseline'));
epochAnaWin.cue = epochAnaParams.anaWins(strcmpi({epochAnaParams.anaWins.name}, 'cue'));
epochAnaWin.delay = epochAnaParams.anaWins(strcmpi({epochAnaParams.anaWins.name}, 'delay'));
epochAnaWin.resp = epochAnaParams.anaWins(strcmpi({epochAnaParams.anaWins.name}, 'response'));

% Feature reduction parameters (using datahigh)
dhParams = struct(...
    'binWidth', analysisParams.binWidth,...
    'use_sqrt', true,...
    'kern', analysisParams.kernSD,...
    'trial_average', false,...
    'keep_neurons', true(1,1));
featureNames = {'gpfa', 'canon_cue', 'canon_delay', 'canon_resp'};
nFeatureSets = length(featureNames);
maxTrajDims = 8;
dist_out = nan(nSessions, nFeatureSets-1);

remapFNames = {...
    'canon_cue' 'canon.,cue';
    'canon_delay' 'canon.,delay';
    'canon_resp' 'canon.,resp.'};
titleNames = featureNames;
for f_ix = 1:length(featureNames)
    if any(strcmpi(featureNames{f_ix}, remapFNames(:,1)))
        r_ix = strcmpi(featureNames{f_ix}, remapFNames(:,1));
        titleNames{f_ix} = remapFNames{r_ix, 2};
    end
end
clear remapFNames f_ix r_ix


sess_ix = 10;  %Best plotting
%%
for sess_ix = 1:nSessions
    %%
    this_sess = sessions(sess_ix);
    thisParams = analysisParams;
    
    %% Load and process ptbmat
    path = 'G:\Projects\MonkeyPFCSaccadeStudies\StudyCentreOut\Data\Preprocessed\ptb';
    data_path = append(path, '\',this_sess,'.ptbmat');
    ptb = load(fullfile(data_path), '-mat');
%     ptb = load(fullfile(paths.preprocessed, 'ptb', this_sess.ptbfname), '-mat');
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
    
    %% Identify neurons with 0 variance across trials in any epoch as 'bad'
    % Do for each CV
    nEpochs = length(epochAnaParams.anaWins);
    nNeurons = length(cDat.invalidUnits);
    goodEpochBool = false(nNeurons, nEpochs);
    for ep_ix = 1:nEpochs
        ana_win = epochAnaParams.anaWins(ep_ix);
        raster = getEpochRaster(cDat, ana_win);  % trials x units x samples
        nGoodSamps = sum(~isnan(raster), 3);
        nSpikes = sum(raster > 0, 3);
        frate = 1000*nSpikes./nGoodSamps;
        goodEpochBool(:, ep_ix) = var(frate) > 0;
    end
    keepBool = ~any(~goodEpochBool, 2);
    clear goodEpochBool ep_ix ana_win raster nGoodSamps nSpikes frate cv_ix
    
    %% Trim trials
    cDat = trimTrials(cDat, trajAnaWin);
    D = cDat2DataHigh(cDat.trial);
    
    %Remove neurons with one or more bad epochs
    for tr_ix = 1:length(D)
        D(tr_ix).data = D(tr_ix).data(keepBool, :);
    end
    
    tVec = trajAnaWin.winEdges(1):trajAnaWin.winEdges(2);
    tVec = tVec(1:end-mod(length(tVec), analysisParams.binWidth));
    tVec = reshape(tVec, analysisParams.binWidth, []);
    step_edges = [tVec(1,:);tVec(end,:)];
    tVec = mean(tVec);
    nSteps = length(tVec);
    
    %% Map each trial's trajectory-bin to the common steps
    % This will be used to align data during MI calculation
    nSamples = cell2mat(cellfun(@size, {D.data},...
        repmat({2}, 1, length(D)), 'UniformOutput', false));
    tr_first_steps = nan(1, nTrials);
    for tr_ix = 1:nTrials
        tr_first_steps(tr_ix) = find(tVec > cDat.trial(tr_ix).tVec(1), 1, 'first');
    end
    tr_last_steps = tr_first_steps - 1 + floor(nSamples/analysisParams.binWidth);
    
    
    
    %% For each feature set
    %fs_ix = find(strcmpi(featureNames, 'canon,cue'));
    dat_out = cell(1, nFeatureSets);
    clear gpfa_out
    for fs_ix = 1:nFeatureSets
        fsNames = strsplit(featureNames{fs_ix}, '_');
        fsName = fsNames{1};
        
        %% Get the trajectories
        
        % Reset the DataHigh parameters
        dhParams.binWidth = analysisParams.binWidth;
        dhParams.kern = analysisParams.kernSD;
        dhParams.keep_neurons = cDat.invalidUnits(keepBool) == DEF.unitState.valid;
        nDimRed = sum(dhParams.keep_neurons);
        
        dhRedMode = -1;
        %smooth: -1; pca:1; fa:3; lda: 4, gpfa:5
        if strcmpi(fsName, 'gpfa')
            dhRedMode = 5;
            dhParams.kern = 0;
            nDimRed = maxTrajDims;
        elseif strcmpi(fsName, 'canon')
            %keep -1 for smooth. TODO: Compare manova1 result with that
            %from DataHigh for LDA when using only baseline epoch
        end
        
        [trajectories, C, lat, dhRedParams] = reducedims(D, dhRedMode, nDimRed, dhParams);
        
        if strcmpi(fsName, 'canon')
            %Get epoch activity
            raster = getEpochRaster(cDat, epochAnaWin.(fsNames{2}));
            raster = raster(:, keepBool, :);
            
            %convert to frate
            nGoodSamps = sum(~isnan(raster), 3);
            nSpikes = sum(raster > 0, 3);
            frate = sqrt(1000*nSpikes./nGoodSamps);
            
            % Calculate canonical weights -> canonical trajectories
            try
                [~, ~, stats] = manova1(frate, classId);
                W = stats.eigenvec(:, 1:maxTrajDims);
                trajectories = centerAndTransformStructArray(trajectories, W);
            catch ME
                trajectories = nan;
                if strcmpi(ME.identifier, 'stats:manova1:SingularSumSquares2')
                    fprintf('MANOVA1 failed for session %i %s\n', sess_ix, fsNames{2});
                else
                    rethrow(ME)
                end
            end
        end
        
        %% Collect trajectories
        if isstruct(trajectories)
            nTrajDim = size(trajectories(1).data, 1);
            nSteps = length(tVec);
            nDims = size(trajectories(1).data, 1);
            nTrials = length(trajectories);
            allDat = nan(nDims, nSteps, nTrials);
            for tr_ix = 1:nTrials
                allDat(:, 1:size(trajectories(tr_ix).data, 2), tr_ix) =...
                    trajectories(tr_ix).data;
            end
            clear tr_ix
        end
        
        %% Compare trajectories between canon and gpfa
        if strcmpi(fsName, 'gpfa')
            gpfa_out = reshape(allDat, [nDims, nSteps*nTrials]);
        elseif strcmpi(fsName, 'canon') && isstruct(trajectories)
            tmp = reshape(allDat, [nDims, nSteps*nTrials]);
            
            if true
                constrainedX = minDistBySwapAndInvert(gpfa_out, tmp);
                
                testDat = constrainedX'*tmp;
                total_dist = sum(sqrt(sum((testDat - gpfa_out).^2, 1)));
                dist_out(sess_ix, fs_ix-1) = total_dist;
                
                % negate some allDat dimensions so they are closer to GPFA;
                % this makes for nicer plots. But don't swap dims.
                [ii,jj,s] = find(constrainedX);
                allDat(ii(s<0), :, :) = -1*allDat(ii(s<0), :, :);
                
            % If you don't believe that the constrainedX is the best way to
            % swap and negate axes to minimize the distance to GPFA, you
            % can uncomment the below code to do an exhaustive search. This
            % is very slow as there are 40320 permutations * 256 ways to
            % negate axes = 10321920 iterations.
            % Note that the below code is incomplete.
            else
                
                p = flipud(perms(1:8));
                nPerms = size(p, 1);
                n = logical(de2bi(0:255));  % Which columns should be negated.
                nNegs = size(n, 1);
                nst = nSteps*nTrials;
                test_dist = Inf;
                perm_ix = 0;
                test_out = nan(nPerms, nNegs);
                while test_dist >= total_dist && perm_ix < nPerms
                    perm_ix = perm_ix + 1;
                    testDat = tmp(p(perm_ix, :), :);
                    for n_ix = 1:nNegs
                        doNeg = n(n_ix, :);
                        nTestDat = testDat;
                        nTestDat(doNeg, :) = -1*nTestDat(doNeg, :);
                        test_out(perm_ix, n_ix) = sum(sqrt(sum((nTestDat - gpfa_out).^2, 1)));
                    end
                    [test_dist, n_ix] = min(test_out);
                    if test_dist == total_dist
                        fprintf('Found matching distance at perm %i neg %i.\n', perm_ix, n_ix);
                    end
                    if mod(perm_ix, 500) == 0
                        fprintf('%i / %i perms...\n', perm_ix, nPerms);
                    end
                end
                
                
            end
        end
        
        
        %% Get mean + SEM trajectories for each class
        if isstruct(trajectories)
            meanSteDat = nan(nDims, nSteps, nClasses, 2);
            for class_ix = 1:nClasses
                class_bool = classId == uqClasses(class_ix);
                meanSteDat(:, :, class_ix, 1) = nanmean(allDat(:, :, class_bool), 3);
                meanSteDat(:, :, class_ix, 2) = nanstd(allDat(:, :, class_bool), 0, 3)...
                    ./ sqrt(sum(~isnan(allDat(:, :, class_bool)), 3));
            end
            
            dat_out{fs_ix} = meanSteDat;
        else
            dat_out{fs_ix} = nan;
        end
        clear class_ix class_bool meanSteDat
    end  %fs_ix
    
    %% Setup the per-session figure
    desiredRes = 300 / 2.54;  %dpi / cm-p-i
    desiredWidth = DEF.fig_width_cm(end) * desiredRes;
    fsize_ax = 18;
    fsize_lab = 24;
    
    myFig = figure('Name', ['Trajectories ' num2str(sess_ix)], ...
        'Position', [10 10 desiredWidth desiredWidth*2/3]);
    colors = colormap('jet');
    colors = colors(round(linspace(1, size(colors, 1), nClasses)), :);
    
    sp_w = 0.225;
    sp_h = 0.26;
    sp_x = 0.032:0.245:1;
    sp_y = flip(0.07:0.31:0.9);
    
    %% Plot the first 3 trajectory dimensions by class
    this_tvec = [tVec'; flip(tVec'); tVec(1)];
    for fs_ix = 1:nFeatureSets
        if ~isnan(dat_out{fs_ix})
            meanSteDat = dat_out{fs_ix};
            ylims = [min(min(min(meanSteDat(1:3, :, :, 1) - meanSteDat(1:3, :, :, 2))))...
                max(max(max(meanSteDat(1:3, :, :, 1) + meanSteDat(1:3, :, :, 2))))];
%             ylims = [-0.8 0.8];
            for dim_ix = 1:3
                %subplot(3, nFeatureSets, nFeatureSets*(dim_ix-1)+fs_ix)
                subplot('Position', [sp_x(fs_ix) sp_y(dim_ix) sp_w sp_h]);
                set(gca, 'LineWidth', 2, 'FontSize', fsize_ax)
                set(gca, 'Color', 'None')
                box off
                hold on
                hp = [];
                for class_ix = 1:nClasses
                    this_dat = squeeze(meanSteDat(dim_ix, :, class_ix, :));
                    this_dat = [this_dat(:, 1)+this_dat(:, 2);...
                        flip(this_dat(:, 1))-flip(this_dat(:, 2));...
                        this_dat(1, 1)+this_dat(1, 2)];
                    hp = [hp, patch(this_tvec, this_dat, colors(class_ix, :))];
                    hp(end).FaceAlpha = 0.5;
                    hp(end).LineStyle = 'none';
                end
                ylim([-max(abs(ylims)) max(abs(ylims))]);
                plot([0 0], get(gca, 'YLim'), 'k--')
                xlim([tVec(1) tVec(end)])
                %set(gca, 'YTick', [])
                hold off
                
                if dim_ix == 1
                    title(strsplit(titleNames{fs_ix},','),...
                        'interpreter', 'tex');
                end
                if dim_ix < 3
                    %set(gca, 'XTick', [])
                    set(gca, 'XTickLabel', {})
                else
                    xlabel('Time After Cue Onset (ms)')
                end
                if fs_ix == 2 && dim_ix == 2
                    hl = legend(hp, classStr,...
                        'Location', 'North',...
                        'Orientation', 'horizontal');
                    hl.Box = 'off';
                    hl.Position = [0.35 0.65 0.3 0.03];
                end
                if fs_ix == 1
                    ylabel(['Dim ' num2str(dim_ix) ' (A.U.)']);
                else
                    set(gca, 'YTickLabel', {});
                end
            end
        end
    end %for each feature set
    
    
    %%
    fprintf(['Finished sessions ' num2str(sess_ix) ' of ' num2str(nSessions) '.\n']);
end

beep

% tbl = table(...
%     reshape(repmat(1:nSessions, 3, 1)', [], 1),...
%     cat(1,...
%         repmat(featureNames(2), nSessions, 1),...
%         repmat(featureNames(3), nSessions, 1),...
%         repmat(featureNames(4), nSessions, 1)), ...
%     dist_out(:),...
%     'VariableNames', {'Session', 'Epoch' 'Distance2GPFA'});

tbl = table((1:nSessions)', dist_out(:, 1), dist_out(:, 2), dist_out(:, 3),...
    'VariableNames', {'Session', 'Epoch1', 'Epoch2', 'Epoch3'});
rm = fitrm(tbl, 'Epoch1-Epoch3 ~ 1');
ranovatbl = ranova(rm);
mc_tbl = multcompare(rm, 'Time', 'ComparisonType', 'bonferroni');
fprintf(...
    ['The average distance between GPFA trajectories and canonical trajectories was ',...
    'significantly different among epochs (F(%i,%i) = %f, p = %f).\n'],...
    ranovatbl{1,2}, ranovatbl{2,2}, ranovatbl{1,4}, ranovatbl{1,5});

[~, ix] = min(nanmean(dist_out));
larg_p = max(mc_tbl.pValue(mc_tbl.Time_1 == ix));
fprintf(['The distance between %s trajectories and GPFA trajectories was ',...
    'less than the distances between the other epoch canonical trajectories ',...
    'and GPFA trajectories (Bonferroni-corrected p <= %f).\n'],...
    featureNames{ix+1}, larg_p);