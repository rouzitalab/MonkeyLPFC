function previewGazePerTrial(this)
%previewGazePerTrial(eData.trial(1))
%Plots the gazeData, the saccades, and some other details.

%% Constants
% According to invalidateGazeData, eye workspace lims in degrees are
degXLim = [-4 44]; degYLim = [-3 33];
aspRat = diff(degXLim)/diff(degYLim);
cmap = colormap('jet');
ncolours = size(cmap, 1);

%TODO: annotate saccades
%http://undocumentedmatlab.com/blog/pinning-annotations-to-graphs/

%% Plot1 - 2D plot.
subplot('position', [aspRat*0.40/2 0.55 aspRat*0.40 0.40])
cla
hold on
axis ij
xlim(degXLim), ylim(degYLim)
titlestr = sprintf('%s Trial %i, Class %i, %i saccade(s).', this.expType, this.expTrial, this.class, this.nSaccades);
title(titlestr);

if ~isempty(this.eyePosition)
    xvec = this.eyePosition(:,1);
    XY = this.eyePosition(:,[2 3]);
    segEdges = [-Inf this.flipScreen Inf];
    nSegs = length(segEdges)-1;
    buf = mod(ncolours, nSegs)/2;
    colInds = round(buf+1:(ncolours/nSegs):ncolours-buf+1);
    
    %Large X's indicate flipScreen events.
    for segIx = 1:nSegs
        xBool = xvec >= segEdges(segIx) & xvec < segEdges(segIx+1);
        if any(xBool)
            thisXY = XY(xBool, :);
            scatter(thisXY(:, 1), thisXY(:, 2), 100, cmap(colInds(segIx), :), '.')
            scatter(thisXY(end, 1), thisXY(end, 2), 300, cmap(colInds(segIx), :), 'x', 'LineWidth', 3)
        end
    end
    
    %Mark identified saccades.
    for sac_ix = 1:length(this.saccades)
        sac = this.saccades(sac_ix);
        plot([sac.startPt(1) sac.endPt(1)], [sac.startPt(2) sac.endPt(2)],...
            'k', 'LineWidth', 2)
        scatter(sac.startPt(1), sac.startPt(2), 300, 'ko')
        scatter(sac.endPt(1), sac.endPt(2), 300, 'kx')
    end
end
hold off

%% Plot2 - X and Y positions over time.
fns = fieldnames(this);
subplot('position', [0.05 0.08 0.9 0.4])
cla
hold on
if ~isempty(this.eyePosition)
    %Plot the eye x and y positions over time.
    xvec = this.eyePosition(:,1); % seconds since the trial started (we probably only took eyePosition after flipScreen(2))
    plot(xvec, XY, 'LineWidth', 1)
    ylims = get(gca, 'YLim');
    
    %Plot vertical lines for every flipscreen.
    vh = [];
    vl = {};
    for segIx = 1:nSegs
        x = segEdges(segIx + 1);
        if x < Inf
            h = plot([x x], ylims, ...
                'LineStyle', '--', ...
                'LineWidth', 1, ...
                'Color', cmap(colInds(segIx), :));
            vh = cat(1, vh, h);
            
            %search for a label that might match.
            thisvl = '';
            for fnIx = 1:length(fns)
                fval = this.(fns{fnIx});
                if ~isempty(fval) && isnumeric(fval) && fval(1) == x
                    thisvl = fns{fnIx};
                end
            end
            if isempty(thisvl) || strcmpi(thisvl, 'flipScreen')
                thisvl = ['flip' num2str(segIx)];
            end
            if length(thisvl) > 4 && strcmpi(thisvl(end-3:end), 'Time')
                thisvl = thisvl(1:end-4);
            end
            vl = cat(1, vl, thisvl);
        end
    end
    legend(vh, vl, 'Location', 'South', 'Orientation', 'horizontal');
    
    %Mark identified saccades.
    for sac_ix = 1:length(this.saccades)
        sac = this.saccades(sac_ix);
        plot(xvec([sac.gazeSampleStart sac.gazeSampleEnd]),...
            [sac.startPt;sac.endPt], 'LineWidth', 3)
        scatter(repmat(xvec(sac.gazeSampleStart),2,1), sac.startPt, 100, 'ko')
        scatter(repmat(xvec(sac.gazeSampleEnd),2,1), sac.endPt, 300, 'kx')
    end
    axis tight
    xlim([min(xvec) max(xvec)])
end
hold off