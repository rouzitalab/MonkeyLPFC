addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths DEF

load(fullfile(paths.results, 'context_ml.mat'));
% 'predictedY', 'actualY', 'actualYstr', 'analysisParams', 'dhParams',...
%     'featureNames', 'mlNames', 'sessions', 'pairBoolOut');

% Because I made a mistake in the original code.
[nSessions, nFeatures] = size(predictedY);
for sess_ix = 1:nSessions
    for feat_ix = 1:nFeatures
        for corner_ix = 1:length(predictedY{sess_ix, feat_ix})
            predictedY{sess_ix, feat_ix}{corner_ix} = ...
                predictedY{sess_ix, feat_ix}{corner_ix}(...
                predictedY{sess_ix, feat_ix}{corner_ix} > 0);
        end
        
    end
end

%% 
[featureNames, xTickNames, re_ix] = normFeatureNames(featureNames);

if any(strcmpi(mlNames, 'SVM-Lin'))
    mlNames{strcmpi(mlNames, 'SVM-Lin')} = 'SVM';
end



%% Calculate accuracy

corner_count = 0;
new_sess_list = [];
for sess_ix = 1:nSessions    
    this_sess_comps = length(predictedY{sess_ix, 1});
    corner_count = corner_count + this_sess_comps;
    new_sess_list = [new_sess_list; repmat(sessions(sess_ix), this_sess_comps, 1)];
end

acc_out = nan(corner_count, nFeatures, length(mlNames));
comp_id = 0;
for sess_ix = 1:size(predictedY, 1)
    for comp_ix = 1:length(actualY{sess_ix});
        trueY = actualY{sess_ix}{comp_ix};
        comp_id = comp_id+1;
        for feat_ix = 1:size(predictedY, 2)
            for ml_ix = 1:length(mlNames)
                predY = predictedY{sess_ix, feat_ix, ml_ix}{comp_ix};
                if size(predY) == size(trueY)
                    acc_out(comp_id, feat_ix, ml_ix) = 100*sum(trueY == predY)/length(predY);
                end
            end
        end
    end
end
acc_out = acc_out(:, re_ix, :);

%%
desiredRes = 300 / 2.54;  %dpi / cm-p-i
desiredWidth = DEF.fig_width_cm(1) * desiredRes;
myfig = figure('Name', 'Classification', 'Units', 'pixels',...
    'Position', [200 200 desiredWidth 2.6*desiredWidth]);
fsize_ax = 18;
fsize_lab = 24;
x_off = 0.1;
f_width = 0.85;
f_height = 0.42;
text_x = -1.0;
text_y = 100;

%% Choose a specific ML, plot and ANOVA2
% Choose rLDA for the comparison between feature sets because feature sets
% with many dimensions don't stand a chance without regularization.
best_ml = find(strcmpi(mlNames, 'SVM'));
target_ML = mlNames{best_ml};
tmp_acc = acc_out(:, :, strcmpi(mlNames, target_ML));

%ranova
sub_block = cell(length(new_sess_list), 1);
for sess_ix = 1:length(new_sess_list)
    sub_block{sess_ix} = [new_sess_list(sess_ix).subject '_' new_sess_list(sess_ix).block];
end
t_id = table(sub_block, 'VariableNames', {'sub_block'});
t_dat = array2table(tmp_acc, 'VariableNames', featureNames);
t = [t_id t_dat];
Meas = table(featureNames','VariableNames',{'Periods'});
rm = fitrm(t,[strjoin(featureNames,',') '~sub_block'],'WithinDesign',Meas);
rtbl = ranova(rm);
tbl = multcompare(rm,'Periods', 'ComparisonType', 'bonferroni')
fprintf(['Baseline (chance) accuracy: %f; ',...
    'repeated measures ANOVA F(%i,%i) = %f, p = %f\n'],...
    mean(tmp_acc(:, strcmpi(featureNames, 'base'))),...
    rtbl{1,2}, rtbl{3,2}, rtbl{1,4}, rtbl{1,5});

%% Plot single-ML acc
%subplot(1,2,1)
symb_size = 200;
subplot('Position', [x_off 0.55 f_width f_height]);
map_bs = struct('A', 'o', 'B', 's', 'C', '^');
map_sf = struct('JerryLee', 'filled', 'Marty', '');
block_symbols = {'o', 's', '^'};
ploth = nan(length(featureNames), length(new_sess_list));
for fs_ix = 1:length(featureNames)
    for sess_ix = 1:length(new_sess_list)
        s = new_sess_list(sess_ix);
        if strcmpi(map_sf.(s.subject), 'filled')
            ploth(fs_ix, sess_ix) = scatter(fs_ix, tmp_acc(sess_ix, fs_ix),...
                symb_size, 'k', map_bs.(s.block), 'filled');
        else
            ploth(fs_ix, sess_ix) = scatter(fs_ix, tmp_acc(sess_ix, fs_ix),...
                symb_size, 'k', map_bs.(s.block));
        end
        hold on
    end
end
set(gca, 'LineWidth', 2)

avg_acc = nanmean(tmp_acc);
bwidth = 0.25;
for fs_ix = 1:length(featureNames)
    plot([fs_ix - bwidth, fs_ix + bwidth], [avg_acc(fs_ix) avg_acc(fs_ix)],...
        'k', 'LineWidth', 3);
end

box off
set(gca, 'FontSize', fsize_ax)
set(gca, 'Color', 'none')
ylim([0 100])
xlim([0 length(xTickNames)+1])
set(gca,'TickLabelInterpreter', 'tex');
set(gca, 'XTick', 1:length(xTickNames))
set(gca, 'XTickLabel', xTickNames)
ylabel([target_ML ' Classification Accuracy (%)'])

sub_ypos = [22 16];
sub_xpos = 5;
sub_fcolor = {'k', 'none'};
sub_str = {'Block A, B, C (JL)', 'Block D, E, F (M)'};
for sub_ix = 1:2
    for bs = 1:length(block_symbols)
        hl = scatter(sub_xpos + (bs-1)*0.5, sub_ypos(sub_ix), symb_size,...
            block_symbols{bs});
        hl.MarkerEdgeColor = 'k';
        hl.MarkerFaceColor = sub_fcolor{sub_ix};
    end
    text(sub_xpos + 1.5, sub_ypos(sub_ix), sub_str{sub_ix},...
        'FontSize', fsize_ax, 'Interpreter', 'tex');
end
hold off
text(text_x, text_y, 'A', 'FontSize', fsize_lab)
%% Plot all results in grouped bars.
subplot('Position', [x_off 0.08 f_width f_height]);
hb = bar(squeeze(nanmean(acc_out)));
hb(1).FaceColor = [1 1 1];
hb(2).FaceColor = [0 0 0];
hb(3).FaceColor = [0.5 0.5 0.5];
box off
set(gca,'TickLabelInterpreter', 'tex');
% set(gca, 'XTick', 1:length(featureNames))
set(gca, 'XTickLabel', xTickNames)
xlim([0 length(featureNames)+1])
ylabel('Classification Accuracy (%)')
ylim([0 100])
hleg = legend(mlNames, 'Location', 'North', 'Orientation', 'horizontal');
hleg.Position = hleg.Position + [0 -0.03 0 0];
legend boxoff
set(gca, 'LineWidth', 2)
set(gca, 'FontSize', fsize_ax)
set(gca, 'Color', 'none')
text(text_x, text_y, 'B', 'FontSize', fsize_lab)

hold on
x = (1:size(acc_out, 2))';
for ml_ix = 1:size(acc_out, 3)
    this_x = x + (ml_ix - 2)*0.23;
    dat = acc_out(:, :, ml_ix);
    y = nanmean(dat);
    y_err = nanstd(dat)./sqrt(size(acc_out, 1));
    heb = errorbar(this_x, y, y_err, 'k', 'LineStyle', 'none');
end
hold off

xlabel('Feature Set')

%%
set(myfig, 'Color', 'none');
savefig(myfig, fullfile(paths.results, 'Figures', 'complement_ml'));