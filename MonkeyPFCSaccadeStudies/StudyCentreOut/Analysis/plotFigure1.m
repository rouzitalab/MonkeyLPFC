function plotFigure1(ptbtrial, cdattrial)
%% Plot for Figure 1
global DEF
desiredRes = 300 / 2.54;  %dpi / cm-p-i
desiredWidth = DEF.fig_width_cm(2) * desiredRes;
myfig = figure('Name', 'Figure1 - Raw Data Example', 'Units', 'pixels',...
    'Position', [200 200 desiredWidth 0.22*desiredWidth]);
fsize_ax = 28;
fsize_lab = 30;

x_off = 0.04;
pl_width = 0.96;
%% Plot eye data
subplot('position', [x_off 0.59 pl_width 0.38])

% xevents = 1000*(ptbtrial.flipScreen(2:5) - ptbtrial.flipScreen(2));
xvec = 1000*(ptbtrial.eyePosition(:,1) - ptbtrial.flipScreen(2));
xlims = 1000*[-0.5 ptbtrial.flipScreen(6)-ptbtrial.flipScreen(2)];
XY = ptbtrial.eyePosition(:,[2 3]);
plot(xvec, XY, 'k', 'LineWidth', 3)

set(gca, 'color', 'none')
box off
set(gca, 'LineWidth', 2)
set(gca, 'FontSize', fsize_ax)
set(gca, 'XTickLabel', {});

tex_x = -250;
[~, xvec_ind] = min(abs(xvec - tex_x));
text(tex_x, XY(xvec_ind, 1)+2, 'X', 'FontSize', fsize_lab)
text(tex_x, XY(xvec_ind, 2)+2, 'Y', 'FontSize', fsize_lab)

xlim(xlims)
ylim([0 40])
s = sprintf('Gaze Pos. (%c)', char(176));
ylabel(s);

hold on
ylims = get(gca, 'YLim');
xevents = [0 300 1300 1600];
for ev_ix = 1:length(xevents)
    plot(xevents([ev_ix ev_ix]), ylims, 'k--', 'LineWidth', 3)
end

%% Plot spike raster
subplot('position', [x_off 0.15 pl_width 0.37])

raster = cdattrial.raster;
[nSamples, nNeurons] = size(raster);

xvec = (1:nSamples) - cdattrial.flipScreen(2);
xvec = repmat(xvec', 1, nNeurons);
xvec = xvec(:);

raster = bsxfun(@times, raster, 1:nNeurons);
raster = raster(:);

scatter(xvec, raster, 200, 'k.');
xlim(xlims);
xlabel('Time After Target Onset (msec)');
ylabel('Neuron ID')
box off
set(gca, 'color', 'none', 'LineWidth', 2, 'FontSize', fsize_ax);

hold on
ylims = get(gca, 'YLim');
for ev_ix = 1:length(xevents)
        plot(xevents([ev_ix ev_ix]), ylims, 'k--', 'LineWidth', 3)
end