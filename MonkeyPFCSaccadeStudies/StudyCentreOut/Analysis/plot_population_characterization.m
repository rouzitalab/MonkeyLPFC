%plot_population_characterization
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths DEF


%% Mutual Information - Not used in the paper
load(fullfile(paths.results, 'mi_trajectories.mat'));
% 'analysisParams', 'trajMiStatP', 'baseMiStatP', 'featureNames');
[nSessions, nWins, nFeatSets] = size(trajMiStatP);


%%
desiredRes = 300 / 2.54;  %dpi / cm-p-i
desiredWidth = DEF.fig_width_cm(1) * desiredRes;
fsize_ax = 14;
fsize_lab = 18;

panelNames = {'A' 'B'};
panelX = {'target onset', 'saccade onset'};
plotPos = {[0.08 0.58 0.9 0.40],[0.08 0.09 0.9 0.40]};


miFig = figure('Name', 'Population Characterization', ...
    'Position', [50 50 desiredWidth desiredWidth/1.5]);
miCmap = colormap(miFig);
col_ix = round(linspace(1, length(miCmap), nFeatSets));

for win_ix = 1:length(analysisParams.anaWins)
    ana_win = analysisParams.anaWins(win_ix);
    subplot('Position', plotPos{win_ix});
    set(gca, 'color', 'none')
    hold on
    
    for fs_ix = 1:nFeatSets
        
        % If this is the dPCA featureset, we need a common xVec across
        % sessions.
        if strcmpi(featureNames{fs_ix}, 'dPCA')
            xrange = [Inf -Inf];
            for sess_ix = 1:nSessions
                sess_xvec = trajMiStatP{sess_ix, win_ix, fs_ix}(:,1);
                xrange = [nanmin([xrange(1);sess_xvec]) nanmax([xrange(2);sess_xvec])];
            end
            xVec = (xrange(1):analysisParams.binWidth:xrange(2))';
            for sess_ix = 1:nSessions
                newdat = [xVec nan(length(xVec), 2)];
                sess_xvec = trajMiStatP{sess_ix, win_ix, fs_ix}(:,1);
                xvec_offset = find(xVec==sess_xvec(1), 1, 'first');
                newdat(xvec_offset:xvec_offset+length(sess_xvec)-1, 2:3) =...
                    trajMiStatP{sess_ix, win_ix, fs_ix}(:, 2:3);
                trajMiStatP{sess_ix, win_ix, fs_ix} = newdat;
            end
        end
        
        % Collapse data across sessions
        allDat = cat(3, trajMiStatP{:, win_ix, fs_ix});
        allMi = squeeze(allDat(:, 2, :));
        miBar = nanmean(allMi, 2);
        miStd = nanstd(allMi, 0, 2) ./ sqrt(sum(~isnan(allMi), 2));
        miPatch = [miBar+miStd; flip(miBar-miStd); miBar(1)+miStd(1)];
        xVec = allDat(:, 1, 1);
        xVecPatch = [xVec; flip(xVec); xVec(1)];
        xVecPatch = xVecPatch(~isnan(miPatch));
        miPatch = miPatch(~isnan(miPatch));
        hp = patch(xVecPatch, miPatch, miCmap(col_ix(fs_ix),:), 'LineStyle', 'none');
        set(hp, 'FaceAlpha', 0.5);
        set(hp, 'EdgeColor', 'none');
        hp_out(fs_ix) = hp;
    end
    hold off
    
    if win_ix == length(analysisParams.anaWins)
        legend(hp_out, featureNames, 'Orientation', 'horizontal', 'Location', 'North')
        legend boxoff
    end
    
    xlim(ana_win.winEdges);
    ylim([0 3]);

    xlabel(['Time relative to ' panelX{win_ix} ' (ms)'])
    box off
    set(gca, 'FontSize', fsize_ax)
    set(gca, 'LineWidth', 2)
    set(gca, 'Color', 'none');
    hold on

    plot([0 0], get(gca, 'YLim'), 'k--')% Plot vertical bar
    ylabel('Mutual Information');

    xl = get(gca, 'XLim');
    yl = get(gca, 'YLim');
    h = text(min(xl) - diff(xl)/20, max(yl), panelNames{win_ix});
    set(h, 'FontSize', fsize_lab);
end
