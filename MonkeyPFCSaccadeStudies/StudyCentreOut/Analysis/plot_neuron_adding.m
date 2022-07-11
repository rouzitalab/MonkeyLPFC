% Plots Figure 10.
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths DEF

%% Load classification. Get baseline epoch accuracy.
load(fullfile(paths.results, 'classification.mat'));
%     'actualY', 'actualYstr', 'analysisParams', 'dhParams',...
%     'featureNames', 'frNames', 'mlNames', 'predictedY', 'sessions');
baseY = predictedY(:, strcmpi(featureNames, 'baseline'), strcmpi(mlNames, 'LDA'));
fullY = predictedY(:, strcmpi(featureNames, 'fullTrial'), strcmpi(mlNames, 'LDA'));
nSessions = length(baseY);
base_acc = nan(nSessions, 1);
full_acc = nan(nSessions, 1);
for sess_ix = 1:nSessions
    base_acc(sess_ix) = 100*sum(baseY{sess_ix} == actualY{sess_ix})./length(actualY{sess_ix});
    full_acc(sess_ix) = 100*sum(fullY{sess_ix} == actualY{sess_ix})./length(actualY{sess_ix});
end

clear baseY fullY anSessions sess_ix
clear actualY actualYstr analysisParams dhParams
clear featureNames frNames mlNames predictedY sessions

%% Load neuron_adding, get accuracy for each neuron addition
load(fullfile(paths.results, 'neuron_adding.mat'));
% 'neur_order', 'neur_acc', 'rneur_acc', 'analysisParams', 'sessions');

nNeurons = max(cellfun(@length, neur_acc));
common_acc = nan(length(neur_acc), nNeurons);
rcommon_acc = nan(length(neur_acc), nNeurons);
for sess_ix = 1:length(neur_acc)
    common_acc(sess_ix, 1:length(neur_acc{sess_ix})) = neur_acc{sess_ix};
    rcommon_acc(sess_ix, 1:size(rneur_acc{sess_ix}, 2)) = nanmean(rneur_acc{sess_ix});
end
clear sess_ix
maxCommonNeurons = find(any(isnan(common_acc)), 2, 'first')-1;
nanmean(common_acc(:, 1:maxCommonNeurons(2)))

%% Compare accuracy at each step to baseline accuracy
p_out_b = nan(1, nNeurons);
t_out_b = nan(1, nNeurons);
df_out_b = nan(1, nNeurons);
rp_out_b = nan(1, nNeurons);
rt_out_b = nan(1, nNeurons);
rdf_out_b = nan(1, nNeurons);

p_out_p = nan(1, nNeurons);
t_out_p = nan(1, nNeurons);
df_out_p = nan(1, nNeurons);
rp_out_p = nan(1, nNeurons);
rt_out_p = nan(1, nNeurons);
rdf_out_p = nan(1, nNeurons);

for neur_ix = 1:nNeurons
    [~, p_out_b(neur_ix), ~, stats] = ttest(common_acc(:, neur_ix), base_acc, 'tail', 'right');
    t_out_b(neur_ix) = stats.tstat;
    df_out_b(neur_ix) = stats.df;
    
    [~, rp_out_b(neur_ix), ~, stats] = ttest(rcommon_acc(:, neur_ix), base_acc, 'tail', 'right');
    rt_out_b(neur_ix) = stats.tstat;
    rdf_out_b(neur_ix) = stats.df;
    
    [~, p_out_p(neur_ix), ~, stats] = ttest(common_acc(:, neur_ix), full_acc, 'tail', 'left');
    t_out_p(neur_ix) = stats.tstat;
    df_out_p(neur_ix) = stats.df;
    
    [~, rp_out_p(neur_ix), ~, stats] = ttest(rcommon_acc(:, neur_ix), full_acc, 'tail', 'left');
    rt_out_p(neur_ix) = stats.tstat;
    rdf_out_p(neur_ix) = stats.df;
end

sig_diff_ix = find(p_out_b < 0.05/nNeurons, 1, 'first');
fprintf(...
    ['Classification of the full trajectories surpassed baseline epoch ',...
    'classification (with all neurons) with the inclusion of the %i ',...
    'best neuron(s) in the model (paired t(%i) = %f, p = %f).\n'],...
    sig_diff_ix, df_out_b(sig_diff_ix), t_out_b(sig_diff_ix), p_out_b(sig_diff_ix)*nNeurons);

sig_diff_ix = find(rp_out_b < 0.05/nNeurons, 1, 'first');
fprintf(...
    ['Classification of the full trajectories surpassed baseline epoch ',...
    'classification (with all neurons) with the inclusion of %i ',...
    'random neuron(s) in the model (paired t(%i) = %f, p = %f).\n'],...
    sig_diff_ix, rdf_out_b(sig_diff_ix), rt_out_b(sig_diff_ix), rp_out_b(sig_diff_ix)*nNeurons);

sig_diff_ix = find(p_out_p > 0.05/nNeurons, 1, 'first');
fprintf(...
    ['Classification of the full trajectories was not worse than peak ',...
    'classification after the inclusion of the %i ',...
    'best neuron(s) in the model (paired t(%i) = %f, p = %f).\n'],...
    sig_diff_ix, df_out_p(sig_diff_ix), t_out_p(sig_diff_ix), p_out_p(sig_diff_ix)*nNeurons);

sig_diff_ix = find(rp_out_p > 0.05/nNeurons, 1, 'first');
fprintf(...
    ['Classification of the full trajectories was not worse than peak ',...
    'classification with the inclusion of %i ',...
    'random neuron(s) in the model (paired t(%i) = %f, p = %f).\n'],...
    sig_diff_ix, rdf_out_p(sig_diff_ix), rt_out_p(sig_diff_ix), rp_out_p(sig_diff_ix)*nNeurons);

clear neur_ix sig_diff_ix
clear p_out_b t_out_b df_out_b
clear rp_out_b rt_out_b rdf_out_b

%TODO: For random neuron adding, compare each session.
% rp_out_sess = cell(length(neur_acc));
% nsig_diff_ix = nan(1, length(neur_acc));
% for sess_ix = 1:length(rneur_acc)
%     rp_out_sess{sess_ix} = nan(1, size(rneur_acc{sess_ix}, 2));
%     for neur_ix = 1:size(rneur_acc{sess_ix}, 2)
%         [~, rp_out_sess{sess_ix}(neur_ix), ~, stats] = ...
%             ttest(rneur_acc{sess_ix}(:, neur_ix), peak_acc(sess_ix), 'tail', 'left');
%     end
%     nsig_diff_ix(sess_ix) = find(rp_out_sess{sess_ix} > 0.05, 1, 'first');
% end

%TODO: Some statement about plateau?

%% Prepare Figure
desiredRes = 300 / 2.54;  %dpi / cm-p-i
desiredWidth = DEF.fig_width_cm(1) * desiredRes;
myfig = figure('Name', 'Neuron Adding', 'Units', 'pixels',...
    'Position', [10 50 desiredWidth desiredWidth*5/3]);
fsize_ax = 18;
fsize_lab = 24;
plotPos = [...
    0.1 0.55 0.88 0.42;...
    0.1 0.08 0.88 0.42;...
    0.8 0.14 0.03 0.04];

%% Plot greedy neuron-adding classification
subplot('Position', plotPos(1,:));
xvec = 1:size(common_acc, 2);
plot(xvec, common_acc', 'k--', 'LineWidth', 1)

% Plot the average + CIs but only up to the number of nuerons common to all
% sessions.
hold on
bAll = sum(isnan(common_acc)) < 3;
muhat = nan(1, length(xvec(bAll)));
muci = nan(2, length(xvec(bAll)));
for step_ix = 1:length(xvec(bAll))
    sess_bool = ~isnan(common_acc(:, step_ix));
    if sum(sess_bool)>1
        [muhat(step_ix), ~, muci(:, step_ix), ~] = ...
            normfit(common_acc(sess_bool, step_ix), 0.05);
    end
end
% goodstep = ~any(isnan(muci));
h = fill([xvec(bAll) fliplr(xvec(bAll)) xvec(1)],...
    [muci(1, :) fliplr(muci(2, :)) muci(1,1)],...
    [0.5 0.5 0.5]);
set(h, 'FaceAlpha', 0.5);
set(h, 'EdgeAlpha', 0);
plot(xvec(bAll), muhat, 'k', 'LineWidth', 3);

box off
ylim([0 100])
xlim([0 size(common_acc, 2)+1])
%set(gca, 'XTickLabel', {});
ylabel({'Greedy Neuron Adding', 'Classification Accuracy (%)'})
set(gca, 'Color', 'none');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', fsize_ax);
text(-0.1*diff(get(gca, 'XLim')), 100, 'A', 'FontSize', fsize_lab)
hold off

%% Plot random neuron-adding classification
subplot('Position', plotPos(2,:));
xvec = 1:size(rcommon_acc, 2);
plot(xvec, rcommon_acc', 'k--', 'LineWidth', 1)

% Plot the average + CIs but only up to the number of nuerons common to all
% sessions.
hold on
bAll = sum(isnan(rcommon_acc)) < 3;
muhat = nan(1, length(xvec(bAll)));
muci = nan(2, length(xvec(bAll)));
for step_ix = 1:length(xvec(bAll))
    sess_bool = ~isnan(rcommon_acc(:, step_ix));
    if sum(sess_bool)>1
        [muhat(step_ix), ~, muci(:, step_ix), ~] = ...
            normfit(rcommon_acc(sess_bool, step_ix), 0.05);
    end
end
% goodstep = ~any(isnan(muci));
h = fill([xvec(bAll) fliplr(xvec(bAll)) xvec(1)],...
    [muci(1, :) fliplr(muci(2, :)) muci(1,1)],...
    [0.5 0.5 0.5]);
set(h, 'FaceAlpha', 0.5);
set(h, 'EdgeAlpha', 0);
plot(xvec(bAll), muhat, 'k', 'LineWidth', 3);

box off
ylim([0 100])
xlim([0 size(common_acc, 2)+1])
xlabel('# Neurons in Model')
ylabel({'Random Neuron Adding', 'Classification Accuracy (%)'})
set(gca, 'Color', 'none');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', fsize_ax);
text(-0.1*diff(get(gca, 'XLim')), 100, 'B', 'FontSize', fsize_lab)
hold off

%% Load neural characterization to get neuron 'type'
load(fullfile(paths.results, 'tuning_epochs.mat'));
% 'trialBoolOut', 'fullInvalidOut', 'tuningCurve',...
% 'tuningStatAN', 'tuningStatMI', 'tuningStatKW',...
% 'tuningStatPopAN', 'tuningStatPopMI',...
% 'analysisParams', 'sessions'
epoch_names = {analysisParams.anaWins.name};
nEpochs = length(epoch_names);
nSessions = length(tuningCurve);
nUnitsPerSession = cellfun(@size, tuningCurve, repmat({1}, nSessions, 1));
nKeptUnits = sum(nUnitsPerSession);
unitStat = cat(1, tuningStatMI{:});
sigTuning = nan(nKeptUnits, nEpochs);
for ep_ix = 1:nEpochs
    sigTuning(:, ep_ix) = unitStat(:, ep_ix, 2) < 0.01;  % Greater than 99% of null.
end

% Classify neurons: tuned or not tuned - baseline, cue, delay, resp (can be more than one)
tuningLabel = 1+sigTuning*[0;4;2;1]; %from 1:8
tuningLabelNames = {'none', 'resp', 'del', 'del&resp', 'cue', 'cue&resp', 'cue&del', 'all'};

%put it back into a per-session structure
tuningLabel = mat2cell(tuningLabel, nUnitsPerSession, 1);

clear epoch_names nEpochs nSessions nUnitsPerSession
clear nKeptUnits unitStat sigTuning ep_ix
clear trialBoolOut fullInvalidOut tuningCurve
clear tuningStatAN tuningStatMI tuningStatKW
clear tuningStatPopAN tuningStatPopMI
clear analysisParams sessions

%% Step through each neuron adding, accumulate labels per class
nNeurons = max(cellfun(@length, neur_order));
neur_labels = nan(length(neur_order), nNeurons);
for sess_ix = 1:length(neur_order)
    neur_ids = neur_order{sess_ix};
    neur_labels(sess_ix, 1:length(neur_order{sess_ix})) = tuningLabel{sess_ix}(neur_ids);
end
clear sess_ix neur_ids

labelCount = nan(nNeurons, length(tuningLabelNames));
for neur_ix = 1:nNeurons
    labelCount(neur_ix, :) = hist(neur_labels(:, neur_ix), 1:length(tuningLabelNames));
end

cumLabelCount = cumsum(labelCount);
cumLabelProportion = bsxfun(@rdivide, cumLabelCount, cumLabelCount(end, :));
res = cumLabelProportion;

%% Compare the order in which neurons were added
% between an in group and an out group

clear grp_cmp

grp_cmp.names = {'all epochs', 'any other combination of epochs'};
grp_cmp.inds = {[8], [2 3 4 5 6 7]};

grp_cmp(end+1) = grp_cmp(end);
grp_cmp(end).names = {'the cue epoch', 'other epochs'};
grp_cmp(end).inds = {[5 6 7 8], [2 3 4]};

grp_cmp(end+1) = grp_cmp(end);
grp_cmp(end).names = {'the delay epoch', 'other epochs'};
grp_cmp(end).inds = {[3 4 7 8], [2 5 6]};

grp_cmp(end+1) = grp_cmp(end);
grp_cmp(end).names = {'the response epoch', 'other epochs'};
grp_cmp(end).inds = {[2 4 6 8], [3 5 7]};

grp_cmp(end+1) = grp_cmp(end);
grp_cmp(end).names = {'only the delay epoch', 'only another epoch'};
grp_cmp(end).inds = {[3], [2 5]};

grp_cmp(end+1) = grp_cmp(end);
grp_cmp(end).names = {'2 or more epochs', '1 or no epochs'};
grp_cmp(end).inds = {[4 6 7 8],[1 2 3 5]};

for cmp_ix = 1:length(grp_cmp)
    is_pref = double(ismember(neur_labels, grp_cmp(cmp_ix).inds{1}));
    is_npref = double(ismember(neur_labels, grp_cmp(cmp_ix).inds{2}));
    is_pref(isnan(neur_labels)) = nan;
    is_npref(isnan(neur_labels)) = nan;
    [~, pref_ind] = find(is_pref == 1);
    [~, nonpref_ind] = find(is_npref == 1);
    [p, h, ~] = ranksum(pref_ind, nonpref_ind);
    [~, typ_ix] = min([median(pref_ind) median(nonpref_ind)]);
    if h == 1
        fprintf(['Neurons with significant tuning during %s were added to ',...
            'the model before neurons with signficant tuning during %s ',...
            '(Wilcoxon rank sum test on selection order, p = %f).\n'],...
            grp_cmp(cmp_ix).names{typ_ix}, grp_cmp(cmp_ix).names{3-typ_ix}, p);
    else
        fprintf('%s vs %s; p = %f\n', grp_cmp(cmp_ix).names{typ_ix}, grp_cmp(cmp_ix).names{3-typ_ix}, p);
    end
end

% Could also try looking at the cumLabelProportion in the next cell

clear has_delay delay_ind nodelay_ind p h stats typ_ix typ_names
clear grp_cmp cmp_ix is_pref is_npref pref_ind nonpref_ind

%%
% For each session, for each step, get the number of neurons of each type
% divided by the total number of neurons of that type in that session
[nSessions, nSteps] = size(neur_labels);
typeList = 1:8;
nTypes = length(typeList);
pOut = nan(nSessions, nSteps, nTypes);
for sess_ix = 1:nSessions
    temp = neur_labels(sess_ix, :);
    type_count = hist(temp, typeList);
    pOut(sess_ix, 1, type_count > 0) = 0;
    for step_ix = 1:nSteps
        this_type = temp(step_ix);
        if ~isnan(this_type)
            pOut(sess_ix, step_ix, this_type) = ...
                sum(temp(1:step_ix) == this_type) / type_count(this_type);
            if step_ix > 1
                pOut(sess_ix, step_ix, setdiff(typeList, this_type)) = ...
                    pOut(sess_ix, step_ix-1, setdiff(typeList, this_type));
            end
        end
    end
end
res = squeeze(nanmean(pOut, 1));
err = squeeze(nanstd(pOut, 0, 1)) ./ squeeze(sqrt(sum(~isnan(pOut), 1)));
clear sess_ix step_ix nSessions nSteps nTypes temp type_count this_type

%%

% Get labels and colours for neuron tuning plot
plotLabelNames = ...
    {'none', 'cue&resp', 'resp', 'del&resp', 'del', 'cue&del', 'cue', 'all'};
plotLabelColours = [...
    128 128 128;... % no tuning = White
    128 128   0;... % cue + resp
    128   0   0;... % resp only = red
    128   0 128;... % delay + resp
      0   0 128;... % delay only = blue
      0 128 128;... % cue + delay
      0 128   0;... % cue only = green
      0   0   0];   % All tuning = Black

%So colours are in more logical order.
[~, remapping] = ismember(plotLabelNames, tuningLabelNames);  %[tuningLabelNames(locb);plotLabelNames]

res = res(:, remapping);
[nSteps, nEpochs] = size(res);

% See if any of the proportions are different to the rest.
pDiff = nan(size(res));
for step_ix = 1:nSteps
    for ep_ix = 1:nEpochs
        other_ep_inds = setdiff(1:nEpochs, ep_ix);
        [~, pDiff(step_ix, ep_ix)] = ...
            ttest(res(step_ix, other_ep_inds), res(step_ix, ep_ix));
    end
end
sigDiff = pDiff*nSteps*nEpochs < 0.05;

subplot('Position', plotPos(2,:))

hp = plot(1:size(res, 1), res,...
    'LineWidth', 3);
hold on
for hp_ix = 1:length(hp)
    hp(hp_ix).Color = plotLabelColours(hp_ix, :)./255;
    
    plot(find(sigDiff(:, hp_ix)), res(sigDiff(:, hp_ix), hp_ix),...
        '*', 'Color', plotLabelColours(hp_ix, :)./255, 'MarkerSize', 20);
end
ylabel({'Proportion of all neurons with significant', 'directional dependence in epoch(s)'})
hold off

% b = bar(cumLabelCount, 'stacked');
% for b_ix = 1:length(b)
%     b(b_ix).FaceColor = plotLabelColours(b_ix, :)./255;
%     if b(b_ix).FaceColor == [1 1 1]
%         b(b_ix).EdgeColor = 'k';
%     else
%         b(b_ix).EdgeColor = 'none';
%     end
% end
% ylabel('# Neurons Across Sessions')

box off
axis tight
xlim([0 size(common_acc, 2)+1])
xlabel('# Neurons in Model')

set(gca, 'Color', 'none');
set(gca, 'LineWidth', 2);
set(gca, 'FontSize', fsize_ax);
text(-0.1*diff(get(gca, 'XLim')), max(get(gca, 'YLim')), 'B', 'FontSize', fsize_lab)

% Add legend
% for cl_ix = 1:length(plotLabelNames)
%     hb = annotation('rectangle',...
%         plotPos(3,:) + (cl_ix-1)*[0 plotPos(3,4) 0 0],...
%         'LineStyle', 'none',...
%         'FaceColor', plotLabelColours(cl_ix, :)./255);
%     if plotLabelColours(cl_ix, :) == [255 255 255]
%         hb.LineStyle = '-';
%         hb.Color = 'k';
%     end
%     
%     hs = annotation('textbox',...
%         plotPos(3,:)+ [plotPos(3,3)+0.005 (cl_ix-1)*plotPos(3,4)-0.01 0 0],...
%         'BackgroundColor', 'none',...
%         'LineStyle', 'none',...
%         'FitBoxToText', 'off',...
%         'Margin', 0,...
%         'FontSize', fsize_ax,...
%         'String', plotLabelNames{cl_ix});
% end
hl = legend(hp, plotLabelNames, 'Location', 'SouthEast');
legend boxoff
hl.FontSize = fsize_ax;


%%
set(myfig, 'Color', 'none');
savefig(myfig, fullfile(paths.results, 'Figures', 'neuron_adding'));