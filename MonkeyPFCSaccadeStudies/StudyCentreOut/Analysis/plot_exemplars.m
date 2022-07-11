addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths DEF
sessions = my_sessions('CentreOut');
analysisParams = my_anaparams('raster');

% % I used the following code to come up with my exemplars.
% load(fullfile(paths.results, 'tuning_epochs.mat'));
% nUnits = cellfun(@size, tuningStatMI, repmat({1}, size(tuningStatMI)));
% sessId = nan(sum(nUnits), 1);
% for sess_ix = 1:length(nUnits)
%     if sess_ix == 1
%         start_ix = 1;
%     else
%         start_ix = sum(nUnits(1:sess_ix-1))+1;
%     end
%     stop_ix = sum(nUnits(1:sess_ix));
%     sessId(start_ix:stop_ix) = sess_ix;
% end
% %'trialBoolOut', 'fullInvalidOut', 'baselineFR', 'tuningCurve',
% %'tuningStatMI', 'analysisParams', 'sessions');
% tuningParams = analysisParams;
% load(fullfile(paths.results, 'dPCA_dat.mat'));
% %'analysisParams', 'nTrials', 'tVec', 'sess_unit_id',...
% %'firingRatesAverage', 'combinedParams', 'margNames', 'W', 'V', 'whichMarg', 'explVar');
% [~, ix] = sort(explVar.margVar(2,:), 'descend');
% comps = ix(1:3);
% neur_weights = sum(abs(W(:, comps)), 2);
% tuningCat = cat(1, tuningStatMI{:});
% for ep_ix = 1:size(tuningCat, 2)
%     subplot(2,2,ep_ix)
%     scatter(neur_weights, tuningCat(:, ep_ix, 1));
%     %Plot unit contribution to time-location vs MI
%     %Choose top units that are also far from origin.
% end
% %Then use the cursor tool to find an x-val associated with an interesting
% %datapoint.
% xval = 1.244;
% [~, ix] = min(abs(neur_weights - xval));
% tmp = cumsum(sessId == sessId(ix));
% exemplar = [sessId(ix) tmp(ix)]

%% Specify plotting details.
plotP.str = {
    'UL' 'UU' 'UR';...  % Mapping the classStr to plot location.
    'LL'  ''  'RR';...
    'DL' 'DD' 'DR'};    % Flipped because y+ve is down
%Position the subplots.
xoff = 2/100; yoff = 2/100; w = 0.3; h = 0.3;
for tmp_ix = 1:numel(plotP.str)
    [i, j] = ind2sub(size(plotP.str), tmp_ix);
    plotP.pos{i, j} = [(j-1)/3+xoff (3-i)/3+yoff w h];
end
clear xoff yoff w h tmp_ix i j
plotP.psth_scale = 1/5;
plotP.psth_size = 1/3;  % Add 1 to the denominator to get portion of graph for psth.
plotP.exemplars = [2 22; 2 25; 5 25; 6 15; 7 8; 10 13; 11 4];  %4 4; 12 1; 4 15; 2 16
%2.22 has high MI in visual that drops off quickly.
%2.25 has high MI in motor but not elsewhere
%5.25 is good delay only
%6.15 has best sustained MI. Also very good weighting.
%7.8 good delay only.
%7.2 interesting delay only, not directional
%10.13 delay only.
%11.4 mostly delay only, different visual/motor

%% Load the data
for sess_ix = unique(plotP.exemplars(:, 1)')
    this_sess = sessions(sess_ix);
    
    %% Load and preprocess the data
    ptb = load(fullfile(paths.preprocessed, 'ptb', this_sess.ptbfname), '-mat');
    [ptb.eData.trial, ptb.params] = getTrialStimInfo(ptb.eData.trial, ptb.params);
    ptb.eData.trial = getTrialBehavResult(ptb.eData.trial, ptb.params, analysisParams);
    ptb.eData.trial = getNewClass(ptb.eData.trial, analysisParams);
    cDat = load(fullfile(paths.preprocessed, 'cDat', [this_sess.edfname(1:end-3) 'mat']));
    cDat = consolidate_ptb_cdat(ptb, cDat);
    trialBool = triageTrials(cDat, analysisParams);
    cDat.trial = cDat.trial(trialBool); clear trialBool
    cDat = triageUnits(cDat, analysisParams);
    
    %% Easy-access, map classes -> directions
    classId = [cDat.trial.newClass];
    [uqClasses, inClassIx] = unique(classId);
    nClasses = length(uqClasses);
    classStr = {cDat.trial(inClassIx).newClassStr};
    classDir = cat(1, cDat.trial(inClassIx).targPol); classDir = classDir(:, 1);
    nUnits = size(cDat.trial(1).raster, 2);
    maxNTrials = max(hist(classId, uqClasses));  % Max number of trials per direction to set y-axis.
    clear inClassIx
    
    %% For the units
    sess_exemplars = plotP.exemplars(plotP.exemplars(:, 1)==sess_ix, 2)';
    for unit_ix = sess_exemplars(sess_exemplars <= nUnits)    
        
        %% Prepare figure for raster plots and tuning curve.
        desiredRes = 300 / 2.54;  %dpi / cm-p-i
        desiredWidth = DEF.fig_width_cm(end) * desiredRes / 2;
        myfig = figure('Name', ['Session ' num2str(sess_ix) ', Unit ' num2str(unit_ix)],...
            'Units', 'pixels',...
            'Position', [200 200 desiredWidth desiredWidth]);
        fsize_ax = 18;
        fsize_lab = 24;

        %Prepare the plots.
        tmpStr = classStr; tmpStr{end+1} = '';
        clear sph
        for sp_ix = 1:(nClasses + 1)
            sph(sp_ix) = subplot('Position', plotP.pos{strcmpi(plotP.str, tmpStr{sp_ix})});
            box off
            set(gca,'Color', 'none')
            set(gca, 'XColor', 'none')
            set(gca, 'YColor', 'none')
            hold on
        end
        clear sp_ix y x tmpStr
        
        %% Get tuningCurve (firing for each analysis window)
        tuningCurve = nan(length(analysisParams.anaWins), nClasses);
        clear hb
        for win_ix = 1:length(analysisParams.anaWins)
            ana_win = analysisParams.anaWins(win_ix);
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
            
            for cl_ix = 1:nClasses
                %%
                % Trials with this class
                trBool = strcmpi({cDat.trial.newClassStr}, classStr{cl_ix});
                
                %Prepare Raster: scatter(spkTimes,spkTrial)
                tmpRaster = raster(trBool, :);
                [spkTrial, ti] = find(tmpRaster > 0);  % trial index and time index of each spike.
                spkTimes = plotXVec(ti);
                [~, tsamp] = find(~isnan(tmpRaster));  % indices of ~nan samples
                sampTimes = plotXVec(tsamp);  % times of ~nan samples
                tuningCurve(win_ix, cl_ix) = 1000*sum(sum(tmpRaster > 0))...
                    /sum(sum(~isnan(tmpRaster)));%Total spiking rate for baseline normalization and tuning curves.
                
                %Prepare PSTH: bar(xEdges, rtSpkBin)
                nSpkBin = histc(spkTimes, xEdges);  % Bin the spike times according to xEdges
                nSampBin = histc(sampTimes, xEdges);  % number of samples per bin
                rtSpkBin = 1000*nSpkBin./nSampBin;  % spike rate (spikes/sec)
                
                % Plot raster
                plot(sph(cl_ix),...
                    spkTimes', spkTrial,...
                    '.', 'MarkerEdgeColor', ana_win.colour);
                
                % Plot the vertical window delimination bar
                if ~isempty(ana_win.vBar)
                    plot(sph(cl_ix),...
                        [ana_win.vBar ana_win.vBar],...
                        [-maxNTrials*plotP.psth_size maxNTrials], 'k--');
                end
                
%                 baseFR = log10(tuningCurve(1, cl_ix));
%                 rtSpkBin = log10(rtSpkBin)...
%                 * (1/baseFR)...  % in units of baseline frate
%                     * plotP.psth_scale...  % Shrink so max = 1/psth_scale
%                     * (maxNTrials*plotP.psth_size);  % How much of the axis are we allowed.
                
                % Plot the PSTH. Axes will be way off for now.
                hb(win_ix, cl_ix) = bar(sph(cl_ix),...
                    xEdges(1:end-1)+analysisParams.binWidth/2,...
                    -1*(log10(rtSpkBin(1:end-1))),...
                    1);
                hb(win_ix, cl_ix).FaceColor = ana_win.colour;
                
                set(sph(cl_ix), 'YDir', 'reverse')  % Invert the Y axis
                set(sph(cl_ix), 'YLim', [-maxNTrials*plotP.psth_size maxNTrials]);
            end
            
            %% Plot tuning curve. Have to fake polar plot around 0,0.
            if ~strcmpi(ana_win.name, 'bcd')
                polRad = [classDir' classDir(1)];
                polRho = [tuningCurve(win_ix, :) tuningCurve(win_ix, 1)];
                [x, y] = pol2cart(polRad, polRho);
                plot(sph(end), x, y, 'LineWidth', 3, 'Color', ana_win.colour);
                clear polRad polRho x y
                axis tight
                axis ij
            end
        end
        
        % Add markings to tuning curve
        plotLim = max(abs([get(sph(end), 'XLim') get(sph(end), 'YLim')]));
        plotLimDiag = pol2cart(deg2rad(45), plotLim);
        set(sph(end), 'XLim', [-plotLim plotLim]);
        set(sph(end), 'YLim', [-plotLim plotLim]);
        plot(sph(end),...
            [-plotLim plotLim], [0 0], 'k--',...
            [0 0], [-plotLim plotLim], 'k--',...
            [-plotLimDiag plotLimDiag], [-plotLimDiag plotLimDiag], 'k--',...
            [-plotLimDiag plotLimDiag], [plotLimDiag -plotLimDiag], 'k--')
        sph(end).XTick = sph(end).XTick(sph(end).XTick > 0);
        for x_ix = 1:length(sph(end).XTick)
            plot(sph(end),...
                [sph(end).XTick(x_ix) sph(end).XTick(x_ix)],...
                [-plotLim/10 plotLim/10], 'k')
            text(sph(end).XTick(x_ix)-plotLim/20, plotLim/8,...
                num2str(sph(end).XTick(x_ix)),...
                'FontSize', fsize_ax);
        end
        text(plotLim/2, plotLim/4, 'Spikes / sec', 'FontSize', fsize_ax)
        
        % Fix PSTH scales
        largest_bar = min(min(arrayfun(@(x) (min(x.YData)), hb)));
        largest_bar = -log10(ceil(10^(-largest_bar)));
        for win_ix = 1:length(analysisParams.anaWins)
            for cl_ix = 1:nClasses
                hb(win_ix, cl_ix).YData = hb(win_ix, cl_ix).YData...
                    ./ largest_bar * -1 * (maxNTrials*plotP.psth_size);
            end
        end
        fprintf('Unit %i-%i PSTH max = %i\n', sess_ix, unit_ix, round(10^(-largest_bar)));
        
        for sp_ix = 1:length(sph)
            set(sph(sp_ix), 'NextPlot', 'replace');
        end
        clear plotLim plotLimDiag sp_ix
    end
    
end












