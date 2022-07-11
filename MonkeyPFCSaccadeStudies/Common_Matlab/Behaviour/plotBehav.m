function [dprime, bias] = plotBehav(tbl, varargin)
%function plotBehav(behav)
% Plots behavioural results
% behav is a table obained from procBehavBlocks

if nargin > 1 && any(strcmpi(tbl.Properties.VariableNames, 'Acc'))
    behavPlotType = varargin{1};
else
    behavPlotType = 'sensitivity';
end

if nargin > 2 &&...
        (any(strcmpi(tbl.Properties.VariableNames, 'DecoderAcc')) ||...
        any(strcmpi(tbl.Properties.VariableNames, 'modelDPrime')) ||...
        any(strcmpi(tbl.Properties.VariableNames, 'DecoderDPrime')))
    y2PlotType = varargin{2};
else
    y2PlotType = 'bias';
end
    
%% Constants
pmarg = 0.01; %If probability approaches 0 or 1, fix by this amount.
arrRadius = 0.05; %Arrow radius in normalized figure units.
arrPosD = [-3 -3]; %Arrow position rel to block middle in units of d'
fsize = 24;

%% Performance Measures
if ~any(strcmpi(tbl.Properties.VariableNames, 'DPrime'))
    [dprime, bias, accuracy] = performanceMeasures(...
        tbl.TrueA, behave.FalseA,...
        tbl.FalseB, tbl.TrueB, pmarg);
else
    dprime = tbl.DPrime;
    bias = tbl.BiasC;
    accuracy = tbl.Acc;
end

%% Plot performance
yyaxis left
if strcmpi(behavPlotType, 'sensitivity')
    h = plot(tbl.TrialIndex, dprime,...
        'LineWidth', 3);
    ylabel('Sensitivity (d'')')
    axis tight
    ylim([-5 5])
    
    hold on
    plot(tbl.TrialIndex, ones(size(dprime))*1.05,...
        'Color', h.Color,...
        'LineStyle', '--',...
        'LineWidth', 1.5)
    hold off
    
elseif strcmpi(behavPlotType, 'accuracy')
    plot(tbl.TrialIndex, accuracy,...
        'LineWidth', 3)
    ylabel('Accuracy (%)')
    axis tight
    ylim([0 100])
end
set(gca, 'LineWidth', 2, 'FontSize', fsize)
set(gca, 'Color', 'none');
box off
xlabel('Block Last Trial Index', 'LineWidth', 2, 'FontSize', fsize)
ylims = get(gca, 'YLim');

yyaxis right
if strcmpi(y2PlotType, 'bias')
    h = plot(tbl.TrialIndex, bias,...
        'LineWidth', 3);
    h.Color = [h.Color 0.3];
    ylabel('Bias',...
        'LineWidth', 2,...
        'FontSize', fsize)
    
    hold on
    plot(tbl.TrialIndex, zeros(size(bias)),...
        'Color', h.Color,...
        'LineStyle', '--',...
        'LineWidth', 1.5)
    hold off
elseif strcmpi(y2PlotType, 'DecoderAcc')
    plot(tbl.TrialIndex, tbl.DecoderAcc,...
        'LineWidth', 3);
    ylabel('Decode Accuracy (%)', 'LineWidth', 2, 'FontSize', fsize)
    if strcmpi(behavPlotType, 'sensitivity')
        ylim([25 100])  % So 70% approx lines up with the d' 1.05 mark.
    else
        ylim([0 100])
    end
elseif strcmpi(y2PlotType, 'modelDPrime')
    plot(tbl.TrialIndex, tbl.modelDPrime,...
        'LineWidth', 3);
    ylabel('modelDPrime', 'LineWidth', 2, 'FontSize', fsize)
elseif strcmpi(y2PlotType, 'DecoderDPrime')
    plot(tbl.TrialIndex, tbl.DecoderDPrime,...
        'LineWidth', 3);
    ylabel('Decoder Sens. (d'')', 'FontSize', fsize)
    ylim([-5 5])
end
box off

%% Plot block divisions
goodStep = tbl.IsGood == 1;
classes = tbl(:,4);
uqBlocks = unique(tbl.BlockIndex(goodStep));

yyaxis left
hold on
for bb=1:length(uqBlocks)
    blockStart = tbl.TrialIndex(find(tbl.BlockIndex==uqBlocks(bb), 1, 'first'));
    plot([blockStart blockStart], ylims, 'k--', 'LineWidth', 1.5);
end
hold off

%% Plot each block's colour-direction rule.
%Plot an arrow for each target/colour combination.
%On top of a white box.
arr_cols = {'r' 'g' 'b'};
arr_ang = [2*pi/4 1*pi/4 0 7*pi/4 6*pi/4 5*pi/4 4*pi/4 3*pi/4];
%8  1   2
%7      3
%6  5   4

for bb = 1:length(uqBlocks)
    bix = uqBlocks(bb);
    blockBool = tbl.BlockIndex==bix;
    tmpBehav = tbl(blockBool, :);
    
    %The arrow origin will be the middle of the block at d=-2 or -4
    arrPosIx = double(mod(bb,2)==0)+1;  % Alternate y position
    arrOrig = [nanmean(tmpBehav.TrialIndex), arrPosD(arrPosIx)];
    [origxf, origyf] = ds2nfu(gca, arrOrig(1), arrOrig(2));
    
    [~, col_ix] = ismember(tmpBehav.CueColor, arr_cols);
    tmp = [tmpBehav.Target col_ix];
    [c, ~, ~] = unique(tmp, 'rows');
    
    for class_ix = 1:length(c)
        [endx, endy] = pol2cart(arr_ang(c(class_ix,1)), arrRadius);
        %TODO: scale endx/endy by the inverse of the figure aspect ratio.
        endx = 3*endx/4;
        endxf = origxf + endx;
        endyf = origyf + endy;
        annotation('arrow', [origxf endxf], [origyf endyf], ...
            'Color', arr_cols{c(class_ix,2)}, 'LineWidth', 5);
    end
end