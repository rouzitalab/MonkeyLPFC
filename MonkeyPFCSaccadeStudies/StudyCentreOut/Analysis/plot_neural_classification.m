addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths DEF

load(fullfile(paths.results, 'classification.mat'));
%     'actualY', 'actualYstr', 'analysisParams', 'dhParams',...
%     'featureNames', 'frNames', 'mlNames', 'predictedY', 'sessions');


%% re-order predictedY
[featureNames, xTickNames, re_ix] = normFeatureNames(featureNames);
predictedY = predictedY(:, re_ix, :);

reorder_mlnames = {'LDA', 'rLDA', 'SVM-Lin'};  %DBN
re_ix = nan(1, length(reorder_mlnames));
for rx = 1:length(reorder_mlnames)
    rebool = strcmpi(mlNames, reorder_mlnames{rx});
    if any(rebool)
        re_ix(rx) = find(rebool);
    end
end
re_ix(isnan(re_ix)) = [];
predictedY = predictedY(:, :, re_ix);
mlNames = reorder_mlnames;
if any(strcmpi(mlNames, 'SVM-Lin'))
    mlNames{strcmpi(mlNames, 'SVM-Lin')} = 'SVM';
end

%% Calculate accuracy
acc_out = nan(length(sessions), length(featureNames), length(mlNames));
sess_ix = 1; fs_ix = 5; ml_ix = 1;
for sess_ix = 1:length(sessions)
    thisY = actualY{sess_ix};
    
    for fs_ix = 1:length(featureNames)
        for ml_ix = 1:length(mlNames)
            predY = predictedY{sess_ix, fs_ix, ml_ix};
            cfnMat = confusionmat(thisY, predY);

            %plotCfn(cfnMat, -90:45:225);
            
            [pCorr, infPerTrial, H] = cfn2Info(cfnMat);
            acc_out(sess_ix, fs_ix, ml_ix) = 100*pCorr;
        end
    end
end

acc_avg = squeeze(mean(acc_out));
acc_avg = acc_avg(:);
[~, best_ix] = max(acc_avg);
[best_ft, best_ml] = ind2sub([length(featureNames), length(mlNames)], best_ix);
fprintf(['Averaged across all sessions, the best machine-learning analysis '...
    'and feature set combination was %s with %s (%f; range %f - %f).\n'],...
    mlNames{best_ml}, featureNames{best_ft}, acc_avg(best_ix),...
    min(acc_out(:, best_ft, best_ml)), max(acc_out(:, best_ft, best_ml)));

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
text_x = -1.3;
text_y = 100;

%% Choose a specific ML, plot and ANOVA2
% Choose rLDA for the comparison between feature sets because feature sets
% with many dimensions don't stand a chance without regularization.
best_ml = find(strcmpi(mlNames, 'rLDA'));
target_ML = mlNames{best_ml};
tmp_acc = acc_out(:, :, strcmpi(mlNames, target_ML));

%ranova
sub_block = cell(length(sessions), 1);
for sess_ix = 1:length(sessions)
    sub_block{sess_ix} = [sessions(sess_ix).subject '_' sessions(sess_ix).block];
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
ploth = nan(length(featureNames), length(sessions));
for fs_ix = 1:length(featureNames)
    for sess_ix = 1:length(sessions)
        s = sessions(sess_ix);
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
sub_xpos = 10;
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
savefig(myfig, fullfile(paths.results, 'Figures', 'classification'));


%% Get predicted saccade end points results from lasso regression
load(fullfile(paths.results, 'correlation.mat'));
% 'actualYnum', 'analysisParams', 'dhParams',...
%     'featureNames', 'predictedYnum', 'sessions', 'targetClass');
[featureNames, xTickNames, re_ix] = normFeatureNames(featureNames);
predictedYnum = predictedYnum(:, re_ix);

targLocs = [...
    512 74;...
    731 165;...
    822 384;...
    731 603;...
    512 694;...
    293 603;...
    202 384;...
    293 165];
targLocs = bsxfun(@rdivide, bsxfun(@minus, targLocs, [512 384]), [512 384]);
nTargs = size(targLocs, 1);

%% Calculate RSq
[nSessions, nFeats] = size(predictedYnum);
nDims = 2;
rOut = nan(nSessions, nFeats, nDims);
accOut = nan(nSessions, nFeats);
rOutClass = nan(nSessions, nFeats, nDims);

for sess_ix = 1:nSessions
    thisY = actualYnum{sess_ix};
    nTrials = size(thisY, 1);
    
    for feat_ix = 1:nFeats
        thisTestY = predictedYnum{sess_ix, feat_ix};
        % Get R^2
        for dim_ix = 1:nDims
            rOut(sess_ix, feat_ix, dim_ix) = ...
                corr(thisY(:, dim_ix), thisTestY(:, dim_ix)).^2;
        end
        
        % Guess target from predicted end point
        dists = nan(nTrials, nTargs);
        for targ_ix = 1:nTargs
            thisTarg = targLocs(targ_ix, :);
            dists(:, targ_ix) = sum(bsxfun(@minus, thisTestY, thisTarg).^2, 2);
        end
        [~, class_ix] = min(dists, [], 2);
        accOut(sess_ix, feat_ix) = ...
            100*sum(targetClass{sess_ix} == class_ix)/nTrials;
        
        % Guess end point from predicted class, correlate with actual end
        % point.
        predClass = predictedY{sess_ix, feat_ix, best_ml};
        predYFromClass = targLocs(predClass, :);
        for dim_ix = 1:nDims
            rOutClass(sess_ix, feat_ix, dim_ix) = ...
                corr(thisY(:, dim_ix), predYFromClass(:, dim_ix)).^2;
        end
    end
       
end
clear sess_ix feat_ix dim_ix targ_ix dists class_ix thisY thisTestY nTrials
clear predClass predYFromClass

[featureNames, xTickNames, re_ix] = normFeatureNames(featureNames);
accOut = accOut(:, re_ix);
rOut = rOut(:, re_ix, :);
rOutClass = rOutClass(:, re_ix, :);

tmp = squeeze(mean(rOut(:, strcmpi(featureNames, 'AvgFRate'), :)));
fprintf(['The average firing rate feature set accounted for'...
    ' %i and %i pcnt of the variance in saccade end point for x and y coordinates, respectively.\n'],...
    round(100*tmp(1)), round(100*tmp(2)));

%% Setup correlations plot
desiredRes = 300 / 2.54;  %dpi / cm-p-i
desiredWidth = DEF.fig_width_cm(1) * desiredRes;
myfig = figure('Name', 'Correlation', 'Units', 'pixels',...
    'Position', [200 200 desiredWidth 2.6*desiredWidth]);
fsize_ax = 18;
fsize_lab = 24;
x_off = 0.1;
f_width = 0.85;
f_height = 0.42;
text_x = -1.3;
text_y = 100;

%% Plot correlations result
subplot('Position', [x_off 0.55 f_width f_height]);
hb = bar(squeeze(nanmean(rOut)));
ylabel('Predicted vs true saccade end-point (R^2)')
hleg = legend({'X', 'Y'}, 'Location', 'North', 'Orientation', 'horizontal');
hleg.Position = hleg.Position + [0 -0.02 0 0];
legend boxoff
set(gca,'TickLabelInterpreter', 'tex');
set(gca, 'XTickLabel', xTickNames)
hb(1).FaceColor = [1 1 1];
hb(2).FaceColor = [0 0 0];
ylim([0 1]);
box off
set(gca, 'color', 'none')
set(gca, 'FontSize', fsize_ax)
set(gca, 'LineWidth', 2)

text(-1.1, max(get(gca, 'YLim')), 'A', 'FontSize', fsize_lab)

hold on
x = (1:size(rOut, 2))';
for dim_ix = 1:size(rOut, 3)
    this_x = x + (dim_ix-1.5)*0.32;
    dat = rOut(:, :, dim_ix);
    y = nanmean(dat);
    y_err = nanstd(dat)./sqrt(size(rOut, 1));
    heb = errorbar(this_x, y, y_err, 'k', 'LineStyle', 'none');
end
hold off

%% Calculate ANOVA on accuracies: Best ML vs from predicted end point.

best_ml = find(strcmpi(mlNames, 'rLDA'));

sess_id = repmat((1:size(accOut, 1))', 1, size(accOut, 2), 2);
feat_id = repmat(1:size(accOut, 2), size(accOut, 1), 1, 2);
meth_id = cat(3, ones(size(accOut)), 2*ones(size(accOut)));
dat = cat(3, accOut, acc_out(:, :, best_ml));

stats = rm_anova2(dat(:), sess_id(:), meth_id(:), feat_id(:), {'Method', 'Feature'});



% fprintf(['Baseline (chance) accuracy: %f; ',...
%     'repeated measures ANOVA F(%i,%i) = %f, p = %f\n'],...
%     mean(tmp_acc(:, strcmpi(featureNames, 'base'))),...
%     rtbl{1,2}, rtbl{3,2}, rtbl{1,4}, rtbl{1,5});


%% Plot accuracies using best, and using predicted end point
subplot('Position', [x_off 0.08 f_width f_height]);

tmp = [squeeze(nanmean(accOut))' squeeze(nanmean(acc_out(:, :, best_ml)))'];
hb = bar(tmp);
ylabel('Target-location prediction accuracy (%)');
hb(1).FaceColor = [1 1 1];
hb(2).FaceColor = [0 0 0];
ylim([0 100]);
box off
set(gca, 'color', 'none')
set(gca, 'FontSize', fsize_ax)
set(gca, 'LineWidth', 2)
hleg = legend({'Using predicted end-point', ['Using ', mlNames{best_ml}]},...
    'Location', 'North', 'Orientation', 'horizontal');
hleg.Position = hleg.Position + [0 -0.02 0 0];
legend boxoff
set(gca, 'TickLabelInterpreter', 'tex');
set(gca, 'XTickLabel', xTickNames);

xlabel('Feature Set')
text(-1.1, max(get(gca, 'YLim')), 'B', 'FontSize', fsize_lab)

hold on
x = (1:size(accOut, 2))';

y = nanmean(accOut);
y_err = nanstd(accOut)./sqrt(size(accOut, 1));
errorbar(x-0.16, y, y_err, 'k', 'LineStyle', 'none');

y = nanmean(acc_out(:, :, best_ml));
y_err = nanstd(acc_out(:, :, best_ml)) ./ sqrt(size(acc_out, 1));
errorbar(x+0.16, y, y_err, 'k', 'LineStyle', 'none');


hold off

%%
set(myfig, 'Color', 'none');
savefig(myfig, fullfile(paths.results, 'Figures', 'correlation'));


%% Plot confusion matrices for each session x featureSet x ML algorithm
%100*sum(bsxfun(@minus, predClass{sess_ix}{fs_ix}, trueClass{sess_ix})==0)./length(trueClass{sess_ix})

% accOutput = nan(size(featureSetNames, 2), length(mlParams), length(sessions));
% for sess_ix = 1:length(sessions)
%     this_sess = sessions(sess_ix);
%     figure('Name', [this_sess.subject ' ' this_sess.block ' ' num2str(this_sess.ntrials)]);
%     nFeatureSets = length(predClass{sess_ix});
%     for fs_ix = 1:nFeatureSets
%         [nTrials, nWins, nML] = size(predClass{sess_ix}{fs_ix});
%         for ml_ix = 1:nML
%             thisPredY = squeeze(predClass{sess_ix}{fs_ix}(:, end, ml_ix));
%             cfnMat = confusionmat(trueClass{sess_ix}, thisPredY);
%             [pCorr, infPerTrial, H] = cfn2Info(cfnMat);
%             plotCfn(cfnMat, -90:45:225);
%             subplot(nFeatureSets, nML, (fs_ix-1)*nML + ml_ix)
%             title([mlParams{ml_ix} ' ' featureSetNames{fs_ix} ' Acc: ' num2str(100*pCorr)]);
%             imagesc(bsxfun(@rdivide, cfnMat, sum(cfnMat, 2)));
%             caxis([0 1]);
%             accOutput(fs_ix, ml_ix, sess_ix) = 100*pCorr;
%         end
%     end
% end
% clear sess_ix this_sess nTrials nFeatureSets nML fs_ix ml_ix
% clear cfnMat pCorr infPerTrial H
% 
% %%
% [nFeatureSets, nML, nSessions] = size(accOutput);
% figure;
% bar3(nanmean(accOutput, 3)');
% set(gca, 'FontSize', 14)
% set(gca, 'YTick', 1:nML);
% set(gca, 'YTickLabel', mlParams);
% ylabel('ML Alg.')
% set(gca, 'XTick', 1:nFeatureSets);
% set(gca, 'XTickLabel', featureSetNames);
% xlabel('Feature Sets')
% zlabel('Accuracy (%)')
% hold on
% h = patch([0 nFeatureSets+1 nFeatureSets+1 0], [0 0 nML+1 nML+1], (100/8)*ones(1,4), 'r');
% set(h, 'FaceAlpha', 0.4);
% hold off
% [p, table, stats] = anova2(reshape(permute(accOutput, [3 1 2]),...
%     nSessions*nFeatureSets, nML), nSessions);


% %% Plot evolution of accuracy over time.
% addpath(fullfile(paths.ml3rd, 'shadedErrorBar'));
% figure;
% cmap = colormap('lines');
% box off
% set(gca, 'LineWidth', 2)
% set(gca, 'FontSize', 20)
% hold on
% for fs_ix = 1:length(featExtractParams)
%     fpara = featExtractParams(fs_ix);
%     xVec = fpara.winEdges(:,2)';
%     accOutput = nan(length(predClass), length(xVec));
%     for sess_ix = 1:length(predClass)
%         predY = predClass{sess_ix}{fs_ix};
%         accOutput(sess_ix, :) = 100*...
%             sum(bsxfun(@minus, predY, trueClass{sess_ix}) == 0)./length(trueClass{sess_ix});
%     end
%     ste = std(accOutput)./sqrt(size(accOutput, 1));
%     lCol = cmap(fs_ix, :);
%     eLine = [mean(accOutput) - ste; mean(accOutput) + ste];
%     xPatch = [xVec, flip(xVec, 2)];
%     yPatch = [mean(accOutput) - ste, flip(mean(accOutput) + ste, 2)];
%     xPatch(isnan(yPatch)) = [];
%     yPatch(isnan(yPatch)) = [];
%     h = patch(xPatch, yPatch, cmap(fs_ix, :), 'edgecolor', 'none', 'facealpha', 0.15);
%     hmain(fs_ix) = plot(xVec, mean(accOutput), 'color', cmap(fs_ix, :), 'LineWidth', 2);
% end
% xlim([min(xVec) max(xVec)])
% xlabel('Data Length (msec)')
% ylabel('SVM Classification Accuracy (%)')
% ylim([0 80])
% legend(hmain, {featExtractParams.fsName}, 'Location', 'NorthWest')
% legend BOXOFF
% plot(get(gca, 'XLim'), [12.5 12.5], 'k--')
% clear fs_ix fpara xVec sess_ix predY ste h hmain  %accOutput

