function plotTuningCurves(seq, varargin)

global paths;
addpath(fullfile(paths.ml3rd, 'mmpolar'));

params.winEdges = [-100 0];
params = varg2params(varargin, params, {'winEdges', 'binWidth', 'kernSD', 'doSqrt', 'doSmooth'});

%% Get each trial's features (i.e. windowed firing)
nTrials = length(seq);
nUnits = size(seq(1).y, 1);
thisSpikes = nan(nTrials, nUnits);
for tr_ix = 1:nTrials
    xBool = seq(tr_ix).tVec >= params.winEdges(1) & seq(tr_ix).tVec < params.winEdges(2);
    thisSpikes(tr_ix, :) = mean(seq(tr_ix).y(:, xBool), 2);
end
% thisSpikes = 100*bsxfun(@rdivide, bsxfun(@minus, thisSpikes, mean(thisSpikes)), mean(thisSpikes));

%% Get each trial's saccade direction.
%     saccades = [cDat.trial(exp_bool).saccades];
%     saccades = [cDat.trial.saccades];
%     sacDir = [saccades.vectDir];
%     sacDir(sacDir<-45/2) = sacDir(sacDir<-45/2)+360;
sacTheta = [seq.sacTheta];
while any(sacTheta < 0)
    sacTheta(sacTheta < 0) = sacTheta(sacTheta < 0) + 2*pi;
end
    
%% Average with bin by sacTheta
targetTheta = deg2rad(0:45:315)';
targetEdges = [targetTheta-(pi/16) targetTheta+(pi/16)];
nTargets = length(targetTheta);
avgCount = nan(nTargets, nUnits);
nPerBin = nan(nTargets, 1);
trialTheta = nan(nTrials, 1);

for bin_ix = 1:size(targetEdges, 1)
    if targetEdges(bin_ix, 1) > targetEdges(bin_ix, 2)
        trialBool = ~(sacTheta > targetEdges(bin_ix, 2) & sacTheta < targetEdges(bin_ix, 1));
    else
        trialBool = sacTheta > targetEdges(bin_ix, 1) & sacTheta < targetEdges(bin_ix, 2);
    end
    trialTheta(trialBool) = bin_ix;
    nPerBin(bin_ix) = sum(trialBool);
    avgCount(bin_ix, :) = mean(thisSpikes(trialBool, :));
end

%% Check for tuning
keepTrialBool = ~isnan(trialTheta);
group = trialTheta(keepTrialBool);
chisq = nan(nUnits, 1);
pvals = nan(nUnits, 1);
pd = nan(nUnits, 1);
for un_ix = 1:nUnits
    x = thisSpikes(keepTrialBool, un_ix);
    [pvals(un_ix), tbl, stats] = kruskalwallis(x, group, 'off');
    chisq(un_ix) = tbl{2,5};
    [~, pd(un_ix)] = max(stats.meanranks);
end

%Print to screen the number of tuned units total, and the number per
%direction
uqgroups = unique(group);
nTuned = nan(length(uqgroups), 1);
for gg_ix = 1:length(uqgroups)
    nTuned(gg_ix) = sum(pd == uqgroups(gg_ix) & pvals < 0.01);
end
fprintf('Found %i tuned units: ', sum(nTuned));
fprintf('%i:%i  ', [rad2deg(targetTheta(uqgroups)) nTuned]');
fprintf('\n');

%% Plot
% Choose up to best 16 units based on kruskall-wallis
% - per-trial frate -> markers
% - avg frate -> red line
figure;
[~, uidx] = sort(pvals, 'ascend');
for uu = 1:min([16 length(uidx)])
    if mod(uu,4) == 0
        xp = 0.75;
    else
        xp = 0.25*(mod(uu,4)-1);
    end
    yp = 1 - 0.25*ceil(uu/4);
    subplot('Position', [xp yp 0.25 0.25])
    uix = uidx(uu);
    
    Y = thisSpikes(:, uix);
    mmpolar(sacTheta', Y', '.', 'Axis', 'off')
    hold on
    
    Y = [avgCount(:, uix); avgCount(1, uix)];
    X = [targetTheta; targetTheta(1)];
    h = mmpolar(X, Y, 'r', 'Axis', 'off', 'TGridVisible', 'on');
    set(h(end), 'LineWidth', 3);
    
%     mmpolar(X, zeros(size(X)), 'k', 'Axis', 'off');
    hold off
end

% %% Fit each unit's firing vs direction to a cosine
% %scatter(sacDir, thisSpikes(:, 1))
% % A*sin(wt + phi) = a*sin(wt) + b*cos(wt) and phi = atan2(a, b)
% X = [sin(sacTheta) cos(sacTheta) ones(size(sacTheta))];
% 
% nUnits = size(thisSpikes, 2);
% amp = nan(1, nUnits);
% phi = nan(1, nUnits);
% bias = nan(1, nUnits);
% stats_out = nan(4, nUnits);
% for uu = 1:nUnits
%     Y = thisSpikes(:, uu);
%     [b, bint, r, rint, stats] = regress(Y, X);  % stats is rsq, f, p
%     bias(uu) = b(3);
%     phi(uu) = rad2deg(atan2(b(2), b(1)));
%     amp(uu) = norm(b(1:2), 2);
%     stats_out(:, uu) = stats;
% end
% while any(phi < 0)
%     phi(phi < 0) = phi(phi < 0) + 360;
% end