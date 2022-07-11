addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths DEF

load(fullfile(paths.results, 'ITR.mat'));
%'actualY', 'actualYstr', 'analysisParams', 'dhParams',...
%    'predictedY', 'sessions', 'winStops');

%% Session accuracy
sess_acc = nan(length(sessions), length(winStops));
for sess_ix = 1:length(sessions)
    sess_acc(sess_ix, :) = 100.*sum(bsxfun(@eq, predictedY{sess_ix}, actualY{sess_ix}))./length(actualY{sess_ix});
end
sess_err = std(sess_acc)./sqrt(length(sessions));
[muhat, sigmahat, muci, sigmaci] = normfit(sess_acc, 0.05);

%% ITR
% Our minimum ITI was about 3 seconds. This includes 700 msec for
% pre-target fixation, 300 for target, 1000 for irrelevant cue, 300 after
% irrelevant cue until fix disappear, then time for saccade reward and
% blank which equals around 700 msec.
% 250 of the 700 pre-target fixation is included in our analysis window. So
% we need to add 450 msec here. We need the 700 msec at the end of the
% trial. So add 1150 to our analysis time.
dur_offset = 1150;
data_dur = dur_offset + winStops - (winStops(1) - dhParams.binWidth);
data_dur = data_dur / 60000; %minutes

%Pierce 1980
%Bits = log2(N_targs) + P_hit * log2(P_hit) + (1-P_hit)*log2[ (1-P_hit) / (N_targs - 1)];
N_targs = 8;
P_hit = sess_acc / 100;
bits = log2(N_targs) + P_hit.*log2(P_hit) + (1-P_hit).*log2( (1-P_hit) / (N_targs - 1) );

bitrate = bsxfun(@rdivide, bits, data_dur);
[brhat, ~, brci, ~] = normfit(bitrate, 0.05);

%% Plot Accuracy

desiredRes = 300 / 2.54;  %dpi / cm-p-i
desiredWidth = DEF.fig_width_cm(1) * desiredRes;
myFig = figure('Name', 'ITR',...
    'Position', [200 200 desiredWidth 2.5*desiredWidth]);
fsize_ax = 18;
fsize_lab = 24;

%%
%subplot(2,1,1)
subplot('Position', [0.11 0.55 0.88 0.43]);
set(gca, 'LineWidth', 2)
plot(winStops, sess_acc', 'color', [0.5 0.5 0.5], 'LineWidth', 0.1)
hold on
h = fill([winStops fliplr(winStops) winStops(1)],...
    [muci(1,:) fliplr(muci(2,:)) muci(1)], [0.5 0.5 0.5]);
set(h, 'FaceAlpha', 0.5);
set(h, 'EdgeAlpha', 0);
plot(winStops, muhat, 'k', 'LineWidth', 3);
hold off
xlims = [min(winStops) - dhParams.binWidth max(winStops)];
xlim(xlims);
set(gca, 'LineWidth', 2)
set(gca, 'Color', 'none');
set(gca, 'XTickLabel', {});
ylim([0 100]);
ylabel('Classification Accuracy (%)');
box off
set(gca, 'FontSize', fsize_ax)
text(min(xlims) - 0.08*abs(diff(xlims)), 100, 'A', 'FontSize', fsize_lab)

%% Plot ITR
%subplot(2,1,2)
symb_size = 200;
map_bs = struct('A', 'o', 'B', 's', 'C', '^');
map_sf = struct('JerryLee', 'filled', 'Marty', '');
block_symbols = {'o', 's', '^'};

axes('position', [0.11 0.09 0.88 0.43]);
plot(winStops, bitrate', 'color', [0.5 0.5 0.5], 'LineWidth', 0.1)
hold on

% Plot *'s at peak bitrate
[~, max_ix] = max(bitrate, [], 2);
ind = sub2ind(size(bitrate), 1:length(sessions), max_ix');
x = winStops(max_ix);
y = bitrate(ind);
for sess_ix = 1:length(sessions)
    s = sessions(sess_ix);
    if strcmpi(map_sf.(s.subject), 'filled')
        scatter(x(sess_ix), y(sess_ix),...
            symb_size, 'k', map_bs.(s.block), 'filled');
    else
        scatter(x(sess_ix), y(sess_ix),...
            symb_size, 'k', map_bs.(s.block));
    end
end

fprintf('\nPeak bitrates:\n');
for sess_ix = 1:length(sessions)
    fprintf('%f bpm @ %i msec\n', bitrate(sess_ix, max_ix(sess_ix)), winStops(max_ix(sess_ix)));
end

% Plot average bitrate
h = fill([winStops fliplr(winStops) winStops(1)],...
    [brci(1,:) fliplr(brci(2,:)) brci(1)], [0.5 0.5 0.5]);
set(h, 'FaceAlpha', 0.5);
set(h, 'EdgeAlpha', 0);
plot(winStops, brhat, 'k', 'LineWidth', 3);
ylims = get(gca, 'YLim');

% Plot | at peak bitrate for average
[~, max_ix] = max(brhat);
plot([winStops(max_ix) winStops(max_ix)], ylims, 'k--', 'LineWidth', 3);
fprintf('\nAveage peak bitrate = %f @ %i\n', brhat(max_ix), winStops(max_ix));

% Plot legend
sub_ypos = [13 8];
sub_xpos = 1500;
sub_fcolor = {'k', 'none'};
sub_str = {'Block A, B, C (JL)', 'Block D, E, F (M)'};
for sub_ix = 1:2
    for bs = 1:length(block_symbols)
        hl = scatter(sub_xpos + (bs-1)*80, sub_ypos(sub_ix), symb_size,...
            block_symbols{bs});
        hl.MarkerEdgeColor = 'k';
        hl.MarkerFaceColor = sub_fcolor{sub_ix};
    end
    text(sub_xpos + 240, sub_ypos(sub_ix), sub_str{sub_ix},...
        'FontSize', fsize_ax, 'Interpreter', 'tex');
end

hold off
set(gca, 'LineWidth', 2)
set(gca, 'Color', 'none');
xlims = [min(winStops) - dhParams.binWidth max(winStops)];
xlim(xlims);
xlabel('Data Segment End (ms after target onset)');
ylabel('Information Transfer Rate (bits / min)');

box off
set(gca, 'FontSize', fsize_ax)
text(min(xlims) - 0.08*abs(diff(xlims)), max(ylims), 'B', 'FontSize', fsize_lab)

%%
set(myFig, 'Color', 'none');
savefig(myFig, fullfile(paths.results, 'Figures', 'ITR'));

%% Apparently two groups of sessions.
u1 = [21 26 27 23];
u2 = [17 22 21 28 14 13 15 20];
[h, p, ci, stats] = ttest2(u1, u2);

t1 = [431-133, 229-33, 258-99, 266-78];
t2 = [259-40, 234-61, 228-46, 328-79,...
    231-25, 208-31, 184-17, 200-34];
[h, p, ci, stats] = ttest2(t1, t2);