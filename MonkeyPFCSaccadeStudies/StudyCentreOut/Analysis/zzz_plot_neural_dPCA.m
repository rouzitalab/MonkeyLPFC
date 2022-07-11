%plot_dPCA
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths DEF
% addpath(fullfile(paths.ml3rd, 'dPCA', 'matlab'));
load(fullfile(paths.results, 'dPCA_dat.mat'));
%'analysisParams', 'nTrials', 'tVec', 'sess_unit_id',...
%'firingRatesAverage', 'combinedParams', 'margNames', 'W', 'V', 'whichMarg', 'explVar');

%% Post process data

%Subtract per-unit mean firing rate.
unitMean = sum(sum(firingRatesAverage, 3), 2)/(size(firingRatesAverage,3)*size(firingRatesAverage,2));
X = bsxfun(@minus, firingRatesAverage, unitMean);
[nNeurons, nClasses, nTimes] = size(X);
nComps = size(W, 2);
X = reshape(X, [nNeurons nClasses*nTimes]);
Z = W' * X;
Z = reshape(Z, [nComps nClasses nTimes]);

%Extract the interesting components for plotting.
nCompsPerRow = 3;  % How many components for each marginilization.
nCompRows = ceil((length(margNames)*nCompsPerRow)/nCompsPerRow);
width = 1/nCompsPerRow;
height = (1/2) / nCompRows;

comp_ix = nan(length(margNames), nCompsPerRow);  % Component number
for marg_ix = 1:length(margNames)
    [~, tmpCompIx] = sort(explVar.margVar(marg_ix, :), 'descend');
    comp_ix(marg_ix, :) = tmpCompIx(1:nCompsPerRow);
end
clear marg_ix tmpCompIx

%Rearrange so time comes first.
marg_ix = repmat((1:length(margNames))', 1, nCompsPerRow);
time_ix = find(strcmpi(margNames, 'Time'));
re_ix = [time_ix setxor(1:length(margNames), time_ix)];
comp_ix = comp_ix(re_ix, :)';
marg_ix = marg_ix(re_ix, :)';
margNames = margNames(re_ix);
clear time_ix re_ix

%% Prepare Plot

margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

time = tVec{1, 1};
time = [time time(end)+50+abs(tVec{1, 2}(1))+tVec{1,2}];

desiredRes = 300 / 2.54;  %dpi / cm-p-i
desiredWidth = DEF.fig_width_cm(2) * desiredRes;
myFig = figure('Name', 'dPCA', ...
    'Position', [200 200 desiredWidth desiredWidth*2/3]);
fsize_ax = 26;
fsize_lab = 36;

tois = [0 250 1300];
toi_symb = {'>' 'o' 's'};
toi_lab = {'target onset' 'delay onset' 'saccade onset'};

% timeEvents = [0 0; 250 250; 1300 1300];
% 
% dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
%     'explainedVar', explVar, ...
%     'marginalizationNames', margNames, ...
%     'marginalizationColours', margColours, ...
%     'whichMarg', whichMarg,                 ...
%     'time', time,                        ...
%     'timeEvents', timeEvents,               ...
%     'timeMarginalization', 2, ...
%     'legendSubplot', 12);

%% Plot components over time
ymax = max(max(max(abs(Z(comp_ix(:), :, :))))); ymax = ymax*1.05;
for c_ix = 1:numel(comp_ix)
    %subplot(4, nCompsPerRow, c_ix)
    col_ix = mod(c_ix-1, nCompsPerRow) + 1;
    row_ix = 1 + floor((c_ix-1)/nCompsPerRow);
    xpos = 1/10 + (col_ix-1)*0.9*width;
    ypos = 1 - row_ix*height;
    subplot('Position', [xpos, ypos, width-1/20, height-0.06]);
    colors = colormap('jet');
    colors = colors(round(linspace(1, size(colors, 1), nClasses)), :);
    
    vaf = round(10*explVar.margVar(marg_ix(c_ix), comp_ix(c_ix)))/10;
    title(['Component #' num2str(comp_ix(c_ix)) ' [' num2str(vaf) '%]']);
    set(gca, 'LineWidth', 2, 'FontSize', fsize_ax)
    set(gca, 'Color', 'None')
    box off
    hold on
    for class_ix = 1:nClasses
        plot(time, squeeze(Z(comp_ix(c_ix), class_ix, :)), 'color', colors(class_ix, :), 'LineWidth', 2)
    end
    for t_ix = 1:length(tois)
        plot([tois(t_ix) tois(t_ix)],[-ymax ymax],'k--');
        if c_ix == 1
            scatter(tois(t_ix), -ymax, 200, 'k', 'filled', toi_symb{t_ix});
        end
    end
    hold off
    xlim([min(time) max(time)])
    set(gca, 'YLim', [-ymax ymax])
    if c_ix == 4
        xlabel('Time (s)')
    else
        set(gca, 'XTickLabel', []);
    end
    if mod(c_ix, 3) == 1
        ylabel({[margNames{marg_ix(c_ix)} ' Comps'] 'Norm Spike Rate'});
        text(1.3*min(time) - 0.3*max(time), ymax, char(65+floor(c_ix/3)), 'FontSize', fsize_lab);
    else
        set(gca, 'YTickLabel', []);
    end
end

%% Plot Legend
%subplot(4, 3, 10)
subplot('Position', [0.3, 0.1, 0.1, height]);
hold on
axis ij
for t_ix = 1:length(tois)
    scatter(0.8, t_ix, 100, 'k', 'filled', toi_symb{t_ix});
    text(1.5, t_ix, toi_lab{t_ix}, 'FontSize', fsize_ax)
end
set(gca, 'Visible', 'off')
classLabels = {'UU' 'UR' 'RR' 'DR' 'DD' 'DL' 'LL' 'UL'};
for cl_ix = 1:length(classLabels)
    plot([0 1], repmat(length(tois) + cl_ix, 1, 2), 'color', colors(cl_ix,:), 'LineWidth', 3);
    text(1.5, length(tois) +  cl_ix, classLabels{cl_ix}, 'FontSize', fsize_ax);
end
xlim([-10 5])
ylim([1 11])

%% Plot VAF
%subplot(4, 3, 7);
subplot('Position', [1/10, 1/14, 1/5, 1/3]);
bar(explVar.margVar', 'stacked', 'BarWidth', 0.75);
box off
set(gca, 'Color', 'none', 'LineWidth', 2, 'FontSize', fsize_ax)
xlim([0 15])
xlabel('Component #')
ylabel('Component Variance (%)')
caxis([0 3])
text(-6, max(get(gca, 'YLim')), 'C', 'FontSize', fsize_lab)

axes('position', [1/6 1/5 0.1 0.1])
d = explVar.totalMarginalizedVar / explVar.totalVar * 100;
pielabels = cell(1, length(margNames));
for m_ix = 1:length(margNames)
    pieLabels{m_ix} = [margNames{m_ix} ' ' num2str(round(d(m_ix))) '%'];
end
hp = pie(d, ones(size(d)), pieLabels);
caxis([0 3])
set(hp(1:2:end), 'LineStyle', 'none')
set(hp(2:2:end), 'FontSize', fsize_ax)
uistack(hp(2:2:end), 'top');


%% Plot trajectories using one marginilization type
plotComps = comp_ix(:, marg_ix(1,:) == find(strcmpi(margNames, 'Location')));
subplot('Position', [0.55 0.07 0.42 0.4]);
hold on
for cl_ix = 1:nClasses
    plot3(squeeze(Z(plotComps(1), cl_ix, :)), squeeze(Z(plotComps(2), cl_ix, :)), squeeze(Z(plotComps(3), cl_ix, :)),...
        'Color', colors(cl_ix, :), 'LineWidth', 2);
    for tt_ix = 1:length(tois)
        time_bool = abs(time-tois(tt_ix)) == min(abs(time-tois(tt_ix)));
        scatter3(Z(plotComps(1), cl_ix, time_bool), Z(plotComps(2), cl_ix, time_bool), Z(plotComps(3), cl_ix, time_bool),...
            200, colors(cl_ix, :), toi_symb{tt_ix}, 'filled');
    end
end
hold off
set(gca, 'LineWidth', 2)
set(gca, 'FontSize', fsize_ax)
set(gca, 'Color', 'none')
xlabel(['Component #' num2str(plotComps(1))]);
ylabel(['Component #' num2str(plotComps(2))]);
zlabel(['Component #' num2str(plotComps(3))]);
text(min(get(gca, 'XLim'))-0.8, max(get(gca, 'YLim')), 'D', 'FontSize', fsize_lab)

%%
set(myFig, 'Color', 'none');
savefig(myFig, fullfile(paths.results, 'Figures', 'dPCA'));

%% Create movie showing trial progress and dPCA component progress.

%% First load some data for eyePosition
sessions = my_sessions('CentreOut');
ptb = load(fullfile(paths.preprocessed, 'ptb', sessions(2).ptbfname), '-mat');
[ptb.eData.trial, ptb.params] = getTrialStimInfo(ptb.eData.trial, ptb.params);
ptb.eData.trial = getTrialBehavResult(ptb.eData.trial, ptb.params, analysisParams);
ptb.eData.trial = getNewClass(ptb.eData.trial, analysisParams);
tr_ix = 73;
eyeDat = ptb.eData.trial(tr_ix).eyePosition;
eyet_Targ = 1000*(eyeDat(:,1) - ptb.eData.trial(tr_ix).flipScreen(2));
eyet_Sac = 1000*(eyeDat(:,1) - ptb.eData.trial(tr_ix).sacStartTime);
eyeDat = eyeDat(:, 2:3);
degXLim = [-4 44]; degYLim = [-3 33];
aspRat = diff(degXLim)/diff(degYLim);
pSize = norm(ptb.params.subjectScreenSize/10) / norm(ptb.params.subjectScreenResolution);
d2m = ptb.params.subjectScreenDistance/10; %100; %distance to monitor in centimeters
targXY = ptb.eData.trial(tr_ix).sacEndXY;
targXY = atand(pSize*targXY/d2m);  %In Deg from top left.
sacBool = false(1, length(time));
sacBool(length(tVec{1,1})+1:end) = true;
sacTime = tVec{1,2};
sacTime = [sacTime(1)+(tVec{1,1} - tVec{1,1}(end))-50 sacTime];

%%
myFig = figure('Position', [200 200 aspRat*800/2 800]);
T = size(Z, 3);
thisZ = Z(plotComps, :, :);
clear F
F(T) = struct('cdata',[],'colormap',[]);
axLims = [min(min(thisZ, [], 3), [], 2) max(max(thisZ, [], 3), [], 2)];

for t_ix = 1:T
    prev_ix = max(1, t_ix-3);
    
    %TODO: Plot task progression
    subplot('position', [0.02 0.52 0.96 0.48])
    cla
    set(gca, 'Color', [0 0 0]);
    xlim(degXLim)
    ylim(degYLim)
    axis ij
    hold on
    if time(t_ix) < 1300
        scatter(mean(degXLim), mean(degYLim), 100, 'ws', 'filled');
    end
    if time(t_ix) > 0
        scatter(targXY(1), targXY(2), 800, 'ws', 'LineWidth', 4);
    end
    eyeb_Targ = eyet_Targ >= time(prev_ix) & eyet_Targ <= time(t_ix);
    eyeb_Sac = eyet_Sac >= sacTime(prev_ix) & eyet_Sac <= sacTime(t_ix);
    eyeb = eyeb_Targ | eyeb_Sac;
    this_eyeDat = eyeDat(eyeb, :);
    plot(this_eyeDat(:,1), this_eyeDat(:,2), 'Color', [0.5 0.5 0.5], 'LineWidth', 3);
    if any(eyeb)
        scatter(this_eyeDat(end,1), this_eyeDat(end,2), 200, colors(2, :), 'o', 'filled');
    end
    hold off
    set(gca, 'XColor', 'none')
    set(gca, 'YColor', 'none')
    set(gca, 'ZColor', 'none')
    
    subplot('position', [0.02 0.02 0.96 0.48])
    cla
    set(gca, 'Color', 'none')
    grayLines = thisZ(:, :, 1:t_ix);
    colLines = thisZ(:, :, prev_ix:t_ix);
    plot3(squeeze(grayLines(1, :, :))', squeeze(grayLines(2, :, :))', squeeze(grayLines(3, :, :))',...
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.5)
    hold on
    for cl_ix = 1:nClasses
        plot3(squeeze(colLines(1, cl_ix, :)), squeeze(colLines(2, cl_ix, :)), squeeze(colLines(3, cl_ix, :)),...
            'Color', colors(cl_ix, :), 'LineWidth', 3);
        
        if time(t_ix) > 0
            targ_ix = find(time > 0, 1, 'first');
            scatter3(...
                squeeze(thisZ(1, cl_ix, targ_ix)),...
                squeeze(thisZ(2, cl_ix, targ_ix)),...
                squeeze(thisZ(3, cl_ix, targ_ix)),...
                200, colors(cl_ix, :), '>', 'filled');
        end
        if time(t_ix) > 250
            del_ix = find(time > 250, 1, 'first');
            scatter3(...
                squeeze(thisZ(1, cl_ix, del_ix)),...
                squeeze(thisZ(2, cl_ix, del_ix)),...
                squeeze(thisZ(3, cl_ix, del_ix)),...
                200, colors(cl_ix, :), 'o', 'filled');
        end
        if sacTime(t_ix) > 0
            sac_ix = find(sacTime > 0, 1, 'first');
            scatter3(...
                squeeze(thisZ(1, cl_ix, sac_ix)),...
                squeeze(thisZ(2, cl_ix, sac_ix)),...
                squeeze(thisZ(3, cl_ix, sac_ix)),...
                200, colors(cl_ix, :), 's', 'filled');
        end
                
    end
    hold off
    set(gca, 'Color', 'none')
    set(gca, 'XColor', 'none')
    set(gca, 'YColor', 'none')
    set(gca, 'ZColor', 'none')
    xlim(axLims(1,:))
    ylim(axLims(2,:))
    zlim(axLims(3,:))
    
    %Add a rotation to the camera
    camorbit(rad2deg(pi*t_ix/T), rad2deg(pi*t_ix/T));
    drawnow
    F(t_ix) = getframe(myFig);
end
%%
writerObj = VideoWriter(fullfile(paths.results, 'traj_behav.mp4'),'MPEG-4');
writerObj.FrameRate = 5;
open(writerObj);
writeVideo(writerObj, F);
close(writerObj);