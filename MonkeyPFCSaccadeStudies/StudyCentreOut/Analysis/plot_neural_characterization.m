% Plots Figure 4. Also prints out a summary of the data.
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths DEF
load(fullfile(paths.results, 'tuning_epochs.mat'));
% 'trialBoolOut', 'fullInvalidOut', 'tuningFRate', 'tuningCountVar',...
%     'tuningCurve', 'sortQuality', 'meanFRate',...
%     'tuningStatAN', 'tuningStatMI', 'tuningStatKW',...
%     'tuningStatPopAN', 'tuningStatPopMI',...
%     'analysisParams', 'sessions');

useUnitStat = 'tuningStatMI'; %'tuningStatAN' 'tuningStatMI' 'tuningStatKW'
useTsStat = 'MI'; %'F' 'MI'
pthresh = 0.001;

%% Classify neurons: tuned or not tuned - baseline, cue, delay, resp (can be more than one)
epoch_names = {analysisParams.anaWins.name};
epoch_names{strcmpi(epoch_names, 'response')} = 'resp';
nEpochs = length(epoch_names);
nSessions = length(tuningCurve);
nUnitsPerSession = cellfun(@size, tuningFRate, repmat({1}, nSessions, 1));
nKeptUnits = sum(nUnitsPerSession);
nClasses = size(tuningCurve{1}, 3);

unitFRate = cat(1, tuningFRate{:});
unitStat = eval(useUnitStat);
unitStat = cat(1, unitStat{:});

sigTuning = nan(nKeptUnits, nEpochs);
for ep_ix = 1:nEpochs
    sigTuning(:, ep_ix) = unitStat(:, ep_ix, 2) < pthresh;  % Greater than 99% of null.
end
sigTuning(isnan(unitStat(:, :, 1))) = nan;
clear ep_ix

%% Collect summary info
nNeurons = size(sigTuning, 1);
nSessions = length(sessions);
bTuneBool = sigTuning(:, strcmpi(epoch_names, 'baseline')) == 1;
zeroFrBool = any(isnan(sigTuning), 2);
keepBool = ~bTuneBool & ~zeroFrBool;
nVisSig = nansum(sigTuning(keepBool, strcmpi(epoch_names, 'cue')));
nVisUq = nansum(sigTuning(keepBool, strcmpi(epoch_names, 'cue')) > 0 ...
    & ~any(sigTuning(keepBool, ~strcmpi(epoch_names, 'cue')), 2));
nDelSig = nansum(sigTuning(keepBool, strcmpi(epoch_names, 'delay')));
nDelUq = nansum(sigTuning(keepBool, strcmpi(epoch_names, 'delay')) > 0 ...
    & ~any(sigTuning(keepBool, ~strcmpi(epoch_names, 'delay')), 2));
nMotSig = nansum(sigTuning(keepBool, strcmpi(epoch_names, 'resp')));
nMotUq = nansum(sigTuning(keepBool, strcmpi(epoch_names, 'resp')) > 0 ...
    & ~any(sigTuning(keepBool, ~strcmpi(epoch_names, 'resp')), 2));

%% Print summary information.
mystr = ['We recorded from %i neurons across %i sessions.'...
    ' Neurons with significant MI in the baseline period (%i)'...
    ' and neurons with zero variance across trials in at least one epoch (%i)'...
    ' were excluded from further single neuron characterization.\n'...
    ' We observed neurons with statistically significant mutual'...
    ' information (MI) among target locations and firing rates during'...
    ' each epoch. Significant MI was found in %i neurons'...
    ' during the cue epoch (%i were significant in only the cue epoch),'...
    ' %i (%i) in the delay epoch, and %i (%i) in the response epoch.\n'];
fprintf(mystr, nNeurons, nSessions, nansum(bTuneBool), sum(zeroFrBool),...
    nVisSig, nVisUq,...
    nDelSig, nDelUq,...
    nMotSig, nMotUq);
clear mystr nVisSig nVisTot nVisUq
clear nDelSig nDelTot nDelUq nMotSig nMotTot nMotUq
clear nBaseline zeroFr nZeroFr nZeroFrAnyEpoch nNeurons

%% Get per session neuron information
allTrials = nan(nSessions, 1);
keepTrials = nan(nSessions, 1);
allUnits = nan(nSessions, 1);
keepUnits = nan(nSessions, 1);
sessTuned = nan(nSessions, nEpochs);

bl_per_bool = strcmpi({analysisParams.anaWins.name}, 'baseline');
other_per_ix = find(~bl_per_bool);

unitStat = eval(useUnitStat);

for sess_ix = 1:nSessions
    allTrials(sess_ix) = length(trialBoolOut{sess_ix});
    keepTrials(sess_ix) = sum(trialBoolOut{sess_ix});
    
    sigTuning = unitStat{sess_ix}(:, :, 2) < pthresh;
    baseTunedBool = sigTuning(:, strcmpi(epoch_names, 'baseline')) == 1;
    zeroVarBool = any(isnan(unitStat{sess_ix}(:, :, 1)), 2);
    keepBool = ~baseTunedBool & ~zeroVarBool;
    
    allUnits(sess_ix) = size(sigTuning, 1);
    keepUnits(sess_ix) = sum(keepBool);
    
    % Tuned during any epoch
    sessTuned(sess_ix, 1) = sum(any(sigTuning(keepBool, other_per_ix), 2));
    % Tuned during specific epoch
    for ep_ix = 2:nEpochs
        sessTuned(sess_ix, ep_ix) = sum(sigTuning(keepBool, ep_ix));
    end
end
clear sess_ix sessInds tuned blTuned per_ix

%% Print a table with summary statistics per session.
colNames = [{'Subject' 'Block' 'AllTrials' 'ExcludedTrials' 'AllNeurons' 'ExcludedNeurons' 'TunedTotal', 'TunedVis', 'TunedDelay', 'TunedResp'}];
myTable = table({sessions.subject}', {sessions.block}',...
    allTrials, allTrials-keepTrials, allUnits, allUnits-keepUnits, sessTuned(:, 1), sessTuned(:, 2), sessTuned(:, 3), sessTuned(:, 4),...
    'VariableNames', colNames);
myTable = sortrows(myTable, {'Subject', 'Block'})

outColNames = {'Subject' 'Block' 'Trials', 'Units', 'Total', 'Cue', 'Delay', 'Resp'};
trStr = cell(height(myTable), 1);
unStr = cell(height(myTable), 1);
for sess_ix = 1:height(myTable)
    trStr{sess_ix} = [num2str(myTable.AllTrials(sess_ix)) ' (' num2str(myTable.ExcludedTrials(sess_ix)) ')'];
    unStr{sess_ix} = [num2str(myTable.AllNeurons(sess_ix)) ' (' num2str(myTable.ExcludedNeurons(sess_ix)) ')'];
end
outTable = table(myTable.Subject, myTable.Block,...
    trStr, unStr, myTable.TunedTotal, myTable.TunedVis, myTable.TunedDelay, myTable.TunedResp,...
    'VariableNames', outColNames);
tblPath = fullfile(paths.results, 'Table1.csv');
writetable(outTable, tblPath, 'FileType', 'text');

clear allTrials keepTrials allUnits keepUnits sessTuned
clear colNames myTable outColnames trStr unStr sess_ix outTable tblPath

%% Plot locational preference across sessions
%per_names = {analysisParams.anaWins.name};
nEpochs = length(epoch_names);
desiredRes = 300 / 2.54;  %dpi / cm-p-i
desiredWidth = DEF.fig_width_cm(end) * desiredRes;
myfig = figure('Name', 'Single Neuron Tuning', 'Units', 'pixels',...
    'Position', [10 50 desiredWidth desiredWidth/3.6]);
fsize_ax = 18;
fsize_lab = 24;

%% Polar plot of average tuning.
if iscell(tuningCurve)
    tuningCurve = cat(1, tuningCurve{:});
end
subplot('position',[50/1500 0.05 400/1500 0.9])
set(gca, 'FontSize', fsize_ax)
polRad = deg2rad([270 315 0 45 90 135 180 225 270]);
clear hlin

unitStat = eval(useUnitStat);
unitStat = cat(1, unitStat{:});

sigTuning = nan(nKeptUnits, nEpochs);
for ep_ix = 1:nEpochs
    sigTuning(:, ep_ix) = unitStat(:, ep_ix, 2) < pthresh;  % Greater than 99% of null.
end
sigTuning(isnan(unitStat(:, :, 1))) = nan;
clear ep_ix

bTuneBool = sigTuning(:, strcmpi(epoch_names, 'baseline')) == 1;
zeroFrBool = any(isnan(sigTuning), 2);
keepBool = ~bTuneBool & ~zeroFrBool;

for ep_ix = 2:nEpochs
    perEpochBool = keepBool & sigTuning(:, ep_ix)>0;
    modIndexBar = nan(1, nClasses);
    modIndexSte = nan(1, nClasses);
    for cl_ix = 1:nClasses
        tmp = tuningCurve(perEpochBool, ep_ix, cl_ix) ./ tuningCurve(perEpochBool, strcmpi(epoch_names, 'baseline'), cl_ix);
        tmpBool = tmp<Inf & tmp>-Inf & ~isnan(tmp);% & tmp~=0;
        modIndexBar(cl_ix) = mean(tmp(tmpBool));
        modIndexSte(cl_ix) = std(tmp(tmpBool))/sqrt(sum(tmpBool));
    end
    inner = modIndexBar - modIndexSte;
    outer = modIndexBar + modIndexSte;
    inner(inner<0) = 0;
    %outer(outer<0) = 0;
    [xin, yin] = pol2cart(polRad, [inner inner(1)]);
    [xout, yout] = pol2cart(polRad, [outer outer(1)]);
    hp = patch([xout,xin], [yout,yin], analysisParams.anaWins(ep_ix).colour, 'linestyle','none', 'facealpha', 0.2);
    hold on
    [x, y] = pol2cart(polRad, [modIndexBar modIndexBar(1)]);
    hlin(ep_ix-1) = plot(x, y, 'LineWidth', 3);
    hlin(ep_ix-1).Color = analysisParams.anaWins(ep_ix).colour;
end

% Plot baseline circle (==1)
hbase = plot(cos(0:0.1:2*pi),sin(0:0.1:2*pi), 'k', 'LineWidth', 3);

% Format plot
axis tight
axis ij
plotLim = max(abs([get(gca, 'XLim') get(gca, 'YLim')]));
plotLimDiag = pol2cart(deg2rad(45), plotLim);
set(gca, 'XLim', [-plotLim plotLim])
set(gca, 'YLim', [-plotLim plotLim])

% Plot directional bars
plot(gca,...
    [-plotLim plotLim], [0 0], 'k--',...
    [0 0], [-plotLim plotLim], 'k--',...
    [-plotLimDiag plotLimDiag], [-plotLimDiag plotLimDiag], 'k--',...
    [-plotLimDiag plotLimDiag], [plotLimDiag -plotLimDiag], 'k--', 'LineWidth', 0.5)

% Plot circular ticks
for r=2:2:6 %[5 10 15 20 25 30]
    plot(r*cos(0:0.1:2*pi),r*sin(0:0.1:2*pi),'k--')
    if mod(r, 2)==0
        text(r, 0.2, num2str(r), 'FontSize', fsize_ax);
    end
end
text(r-2, 0.8, 'Relative Firing Rate', 'FontSize', fsize_ax);
box off


box off
set(gca,'Color', 'none')
set(gca, 'XColor', 'none')
set(gca, 'YColor', 'none')
legend_entries = epoch_names;
legend_entries(strcmpi(legend_entries, 'baseline')) = {'baseline(==1)'};
hleg = legend([hbase hlin], legend_entries, 'Location', 'NorthEast');
legend boxoff
set(hleg, 'FontSize', fsize_ax, 'LineWidth', 2)
h = text(min(get(gca, 'XLim')), min(get(gca, 'YLim')), 'A');
set(h, 'FontSize', fsize_lab)
clear polRad per_ix perEpochBool modIndexBar modIndexSte cl_ix tmp tmpBool
clear inner outer xin yin xout yout hp x y plotLim plotLimDiag hlin hleg

%% Polar plot, lines of length nTuned, coloured for each epoch.

%Preferred locations for each neuron*epoch
[~, prefLoc] = nanmax(tuningCurve, [], 3);

bTuneBool = sigTuning(:, strcmpi(epoch_names, 'baseline')) == 1;
zeroFrBool = any(isnan(sigTuning), 2);
keepBool = ~bTuneBool & ~zeroFrBool;

%Check to see if prefLoc is sig diff to non-pref targets
prefIsDiff = false(size(prefLoc));
tempMC = cat(1, tuningStatANMC{:});  % 1 location vs rest t-test
[nNeurons, nEpochs, nClasses] = size(tempMC);
for neur_ix = 1:nNeurons
    for ep_ix = 1:nEpochs
        prefIsDiff(neur_ix, ep_ix) = ...
            (abs(tempMC(neur_ix, ep_ix, prefLoc(neur_ix, ep_ix))) < 0.05/nNeurons)...
            && (sigTuning(neur_ix, ep_ix)>0);
    end
end

nNeurEps = sum(sum(sigTuning(keepBool, 2:end), 2));
nPrefNeurEps = sum(sum(prefIsDiff(keepBool, 2:end)));

fprintf(['Of the %i neuron-epochs with significant tuning, %i (%f pcnt) ',...
    'met the criterion that the preferred location must be associated with ',...
    'different firing than all other locations.\n'],...
    nNeurEps, nPrefNeurEps, 100*nPrefNeurEps/nNeurEps);

stat_lvr = nan(nEpochs, 2);
for ep_ix = 1:nEpochs
    thisPrefDir = prefLoc(prefIsDiff(:, ep_ix), ep_ix);
    LvR = nan(size(thisPrefDir));
    LvR(ismember(thisPrefDir, [6 7 8])) = 1;
    LvR(ismember(thisPrefDir, [2 3 4])) = 2;
%     LvR(ismember(thisPrefDir, [1 5])) = 3;
    LvR = LvR(~isnan(LvR));
    [h,p,st] = chi2gof(LvR,'ctrs',[1 2],'expected',size(LvR, 1) * [3/6 3/6]); %[3/8 3/8 2/8]);
    stat_lvr(ep_ix, :) = [st.chi2stat p];
    fprintf([analysisParams.anaWins(ep_ix).name, ': x^2 = ',...
        num2str(stat_lvr(ep_ix, 1)), ', p = ', num2str(stat_lvr(ep_ix, 2)),...
        '. LvR=', num2str(sum(LvR==1)), '-' num2str(sum(LvR==2)), '\n']);
end
clear prefDir stat_lvr per_ix thisPrefDir LvR h p st

%% Plot distribution of preferred directions among significant MI neurons.
epAngSpacing = 12;
subplot(1,3,2)
subplot('position',[500/1500 0.05 400/1500 0.9])
%Spread the different epochs over a range of 15 degrees.
classTargCentres = [270 315 0 45 90 135 180 225];
degPerEpoch = epAngSpacing/(nEpochs-1);
epochDegOffsets = flip(abs(-epAngSpacing:degPerEpoch:-0.1));
epochDegOffsets = epochDegOffsets - mean(epochDegOffsets);
maxNTuned = 0;

for ep_ix = 2:nEpochs
    pd = prefLoc(prefIsDiff(:, ep_ix), ep_ix);
    nTuned = hist(pd, 1:nClasses);
    for cl_ix = 1:nClasses
        [x,y] = pol2cart(deg2rad([1;1]*classTargCentres(cl_ix)+epochDegOffsets(ep_ix-1)), [0;nTuned(cl_ix)]);
        plot(x, y, 'Color', analysisParams.anaWins(ep_ix).colour, 'LineWidth', 3),
        hold on
    end
    maxNTuned = max(maxNTuned, nTuned);
end
axis tight
axis ij
plotLim = max(abs([get(gca, 'XLim') get(gca, 'YLim')]));
set(gca, 'XLim', [-plotLim plotLim])
set(gca, 'YLim', [-plotLim plotLim])

% Plot circular ticks
for r=[5 10 15 20] %[5 10 15 20 25 30]
    plot(r*cos(0:0.1:2*pi),r*sin(0:0.1:2*pi),'k--')
    if mod(r, 10)==0
        text(r, -2, num2str(r), 'FontSize', fsize_ax);
    end
end
text(12, -5, '# Neurons', 'FontSize', fsize_ax);
box off

% Cleanup axes
set(gca,'Color', 'none')
set(gca, 'XColor', 'none')
set(gca, 'YColor', 'none')
set(gca, 'FontSize', fsize_ax)
h = text(min(get(gca, 'XLim')), min(get(gca, 'YLim')), 'B');
set(h, 'FontSize', fsize_lab)
clear degPerEpoch epochDegOffsets maxNTuned per_ix pd nTuned cl_ix x y plotLim r

%% Delta angle in PD across periods
%show each neuron within a column, one column for each pairwise comb.
subplot(1,3,3)
subplot('position',[990/1500 0.12 490/1500 0.83])
classTargCentres = [270 315 0 45 90 135 180 225];
nCombs = nchoosek(nEpochs-1,nEpochs-2);
myLabel = cell(1, nCombs);
nOverlap = nan(1, nCombs);
xd = flip(180:-45:-179);
deltaOut = nan(length(xd), nCombs);
comb_ix = 0;
for ep_ix = 2:nEpochs
    for op_ix = 3:nEpochs
        if ep_ix < op_ix
            comb_ix = comb_ix + 1;
            myLabel{comb_ix} = [epoch_names{ep_ix} ' - ' epoch_names{op_ix}];
            
            bothBool = nansum(sigTuning(:, [ep_ix op_ix]), 2) == 2;
            bothBool = bothBool & sum(prefIsDiff(:, [ep_ix op_ix]), 2)==2;
            
            nOverlap(comb_ix) = sum(bothBool);
            deltaAngles = diff(classTargCentres(prefLoc(bothBool, [ep_ix op_ix])), 1, 2);
            deltaAngles(deltaAngles<=-180) = deltaAngles(deltaAngles<=-180) + 360;
            deltaAngles(deltaAngles>180) = deltaAngles(deltaAngles>180) - 360;
            deltaOut(:, comb_ix) = hist(deltaAngles, xd);
        end
    end
end
bar(xd, deltaOut);
xlim([-155 200]);
set(gca,'Color', 'none')
set(gca, 'FontSize', fsize_ax)
set(gca, 'LineWidth', 2)
box off
xlabel('\Delta Angle (^o)')
ylabel('# Neurons')
hl = legend(myLabel);
set(hl, 'FontSize', fsize_ax)
legend boxoff
h = text(min(get(gca, 'XLim'))-0.1*diff(get(gca, 'XLim')),...
    max(get(gca, 'YLim')), 'C');
set(h, 'FontSize', fsize_lab)

fprintf('%f pcnt of significant MI epoch-pairs had the same PD.\n', 100*sum(deltaOut(xd==0,:))/sum(sum(deltaOut)));

clear nCombs myLabel nOverlap xd deltaOut comb_ix per_ix op_ix bothBool deltaAngles hl h

%% Save figure
set(myfig, 'Color', 'none');
savefig(myfig, fullfile(paths.results, 'Figures', 'tuning_average'));
clear myfig

%% Plot firing rate vs tuned epoch(s)

% % Get unit tuning
% sigTuning = nan(nKeptUnits, nEpochs);
% for ep_ix = 1:nEpochs
%     sigTuning(:, ep_ix) = unitStat(:, ep_ix, 2) < pthresh;  % Greater than 99% of null.
% end
% sigTuning(isnan(unitStat(:, :, 1))) = nan;
% 
% bTuneBool = sigTuning(:, strcmpi(epoch_names, 'baseline')) == 1;
% zeroFrBool = any(isnan(sigTuning), 2);
% keepBool = ~bTuneBool & ~zeroFrBool;
% 
% tuningfr = cat(1, tuningFRate{:});
% tuningfr = tuningfr(keepBool, :);
% 
% squal = cat(1, sortQuality{:});
% squal = squal(keepBool);
% 
% sigTuning = logical(sigTuning(keepBool, :));
% 
% clear keepBool zeroFrBool bTuneBool ep_ix
% 
% % There aren't enough neurons with exclusive tuning in cue epoch to look at
% % the effect of exclusive-tuning-epoch on firing rate
% 
% % Look at effect of non-exclusive tuning epoch on firing rate
% bldat = [...
%     tuningfr(sigTuning(:, 2), 1);...
%     tuningfr(sigTuning(:, 3), 1);...
%     tuningfr(sigTuning(:, 4), 1)];
% dat = [...
%     tuningfr(sigTuning(:, 2), 2);...
%     tuningfr(sigTuning(:, 3), 3);...
%     tuningfr(sigTuning(:, 4), 4)];
% grouping = [
%     1*ones(sum(sigTuning(:, 2)), 1);...
%     2*ones(sum(sigTuning(:, 3)), 1);...
%     3*ones(sum(sigTuning(:, 4)), 1)];
% groupNames = {'Cue', 'Delay', 'Resp.'};
% 
% [p, anovatab, stats] = kruskalwallis(bldat, grouping);
% % No effect of tuning epoch(s) on baseline frate
% 
% [p, anovatab, stats] = kruskalwallis(sqrt(dat), grouping);
% % Strong effect of tuning epoch on firing rate (also modulation) within that epoch
% 
% dat = [...
%     squal(sigTuning(:, 2));...
%     squal(sigTuning(:, 3));...
%     squal(sigTuning(:, 4))];
% [p, anovatab, stats] = kruskalwallis(dat, grouping);
% % No effect of tuned epoch on sort quality
% % suggesting that sort quality does not reveal any information about what
% % part of the task the neuron participates in.
% 
% %
% % Look at firing rate vs number of tuned epochs
% %
% 
% %There are 8 possible tuning outcomes
% tunedNames = {'none' 'cue', 'delay', 'resp', 'cue&delay', 'cue&resp', 'delay&resp', 'all3'};
% tunedCodes = [...
%     0 0 0 0;...
%     0 1 0 0;...
%     0 0 1 0;...
%     0 0 0 1;...
%     0 1 1 0;...
%     0 1 0 1;...
%     0 0 1 1;...
%     0 1 1 1];
% nTuned = [0 1 1 1 2 2 2 3];  % Number of epochs with sig tuning
% fr_scaling = [0 1 4 1];  % delay epoch is 4 times longer than other epochs
% 
% % Tuning in each epoch between singly-, doubly-, and triply-tuned neurons
% [~, tuningGroup] = ismember(sigTuning, tunedCodes, 'rows');
% 
% epoch_goi = [2 5 6 8; 3 5 7 8; 4 6 7 8];
% for ep_ix = 1:size(epoch_goi, 1)
%     this_g = epoch_goi(ep_ix, :);
%     [this_bool, this_tungr] = ismember(tuningGroup, this_g);
%     this_fr = tuningfr(this_bool, ep_ix);
%     
%     this_n = nTuned(this_g);
%     this_n = this_n(this_tungr(this_bool))';
%     
% %     [p, annovatab, stats] = kruskalwallis(sqrt(this_fr), this_n);
%     
% %     [p, annovatab, stats] = kruskalwallis(sqrt(this_fr), squal(this_bool));
%     
%     [p, annovatab, stats] = kruskalwallis(squal(this_bool), this_n);
%     
% end
% 
% %
% % For multiply-tuned neurons, look if there is an effect rotated tuning
% % true/false on sort quality.
% %
% 
% squal = cat(1, sortQuality{:});
% 
% % Only consider neurons that have significant preferred target locations.
% [~, prefLoc] = nanmax(tuningCurve, [], 3);
% 
% %First, check to see if prefLoc is sig diff to non-pref targets
% sigTuning = nan(nKeptUnits, nEpochs);
% for ep_ix = 1:nEpochs
%     sigTuning(:, ep_ix) = unitStat(:, ep_ix, 2) < pthresh;  % Greater than 99% of null.
% end
% sigTuning(isnan(unitStat(:, :, 1))) = nan;
% prefIsDiff = false(size(prefLoc));
% tempMC = cat(1, tuningStatANMC{:});  % 1 location vs rest t-test
% [nNeurons, nEpochs, nClasses] = size(tempMC);
% for neur_ix = 1:nNeurons
%     for ep_ix = 1:nEpochs
%         prefIsDiff(neur_ix, ep_ix) = ...
%             (abs(tempMC(neur_ix, ep_ix, prefLoc(neur_ix, ep_ix))) < 0.05/nNeurons)...
%             && (sigTuning(neur_ix, ep_ix)>0);
%     end
% end
% 
% isRotated = nan(nNeurons, 1);
% for neur_ix = 1:nNeurons
%     this_b = prefIsDiff(neur_ix, :);
%     if any(this_b)
%         this_pl = prefLoc(neur_ix, :);
%         isRotated(neur_ix) = length(unique(this_pl(this_b))) > 1;
%     end
% end
% useBool = ~isnan(isRotated);
% [p, annovatab, stats] = kruskalwallis(squal(useBool), isRotated(useBool));

%% Mutual Information

% Setup the figure
desiredRes = 300 / 2.54;  %dpi / cm-p-i
desiredWidth = DEF.fig_width_cm(end) * desiredRes;
miFig = figure('Name', 'Mutual Information', ...
    'Position', [1 50 desiredWidth desiredWidth/3.4]);
fsize_ax = 18;
fsize_lab = 24;

%Change the epoch names.
analysisEpochs = analysisParams.anaWins;  % From the -period analysis.
analysisEpochs(strcmpi({analysisEpochs.name}, 'cue')).name = 'cue';
analysisEpochs(strcmpi({analysisEpochs.name}, 'response')).name = 'resp';

% Whether we are doing normal MI, relative-MI or delta-MI
doRelMI = false;
doDeltaMI = false;
if doRelMI
    baseMI = unitStat(:, strcmpi({analysisEpochs.name}, 'baseline'), 1);
else
    baseMI = ones(size(unitStat(:, 1, 1)));
end

%% MI on epochs
subplot('Position', [0.05 0.1 0.30 0.8]);
set(gca, 'LineWidth', 2);
per_cols = {'w' 'g' 'b' 'r'};
barh = nan(1, length(analysisEpochs));
for ep_ix = 1:length(analysisEpochs)
    ep_dat = unitStat(:, ep_ix, 1)./baseMI;
    if doDeltaMI
        ep_dat = ep_dat - unitStat(:, 1, 1);
    end
    barh(ep_ix) = bar(ep_ix, nanmean(ep_dat), per_cols{ep_ix});
    hold on
    thisErr = nanstd(ep_dat)./sqrt(nKeptUnits);
    errorbar(ep_ix, nanmean(ep_dat), thisErr, thisErr, 'k', 'LineWidth', 2);
end
%TODO: Fix y-axis size and ticks.
if doRelMI
    ylim([0 5]);
elseif doDeltaMI
    ylim([-0.1 1]);
else
    ylim([0 0.4]);
end
hold off;
box off;
set(gca, 'XTick', 1:length(analysisEpochs), 'XTickLabel', {analysisEpochs.name}, 'FontSize', fsize_ax);
set(gca, 'LineWidth', 2)

if doRelMI
    ylabel('Relative MI')
elseif doDeltaMI
    ylabel('\Delta MI (bits)')
else
    ylabel('Mutual Information (bits)')
end

xlabel('Epoch');
set(gca, 'Color', 'none');
text(min(get(gca, 'XLim')) - diff(get(gca, 'XLim'))/7, max(get(gca, 'YLim')), 'A', 'FontSize', fsize_lab)

%% Compare MI across epochs for all neurons

anaNames = {analysisEpochs.name};
% build ranova multcompare p-value confusion matrix
t = array2table([(1:nKeptUnits)', sqrt(unitStat(:, :, 1))],...
    'VariableNames', ['Neuron' anaNames]);
Meas = table(anaNames','VariableNames',{'Epochs'});
rm = fitrm(t,[strjoin(anaNames,',') '~Neuron'],'WithinDesign',Meas);
rtbl = ranova(rm);
tbl = multcompare(rm,'Epochs', 'ComparisonType', 'bonferroni');
cfnP = nan(length(analysisEpochs));
for ep1 = 1:length(analysisEpochs)
    for ep2 = 1:length(analysisEpochs)
        if ep1 < ep2
            cfnP(ep1, ep2) = tbl.pValue(strcmpi(tbl.Epochs_1, anaNames{ep1}) & strcmpi(tbl.Epochs_2, anaNames{ep2}));
        end
        
        if ep2-ep1 == 1
        elseif ep2-ep1 == 2
            if ep1==1
            elseif ep1==2
            end
        elseif ep2-ep1 ==3
        end
    end
end
num2str(cfnP)

%Print statement
fprintf(['The temporal dependence of MI is shown in Figure 5. ', ...
    'Averaged across neurons, MI was significantly greater ', ...
    'during cue, delay, and response epochs than during the ',...
    'baseline epoch (paired t-tests, Bonferroni-corrected p << 0.001).',...
    'Further, MI was significantly different between each combination ',...
    'of epochs, and the delay epoch was associated with the greatest MI',...
    '(Figure 5 A; all Bonferroni-corrected p-values < 0.01).\n']);

%TODO: Add significance bars comparing different epochs

%% Timeseries

load(fullfile(paths.results, 'mi_timeseries.mat'));
% 'analysisParams',...
%     'timeseriesUnitMiP', 'timeseriesUnitANP', ...
%     'timeseriesUnitAN2', 'timeseriesUnitRMA', ...
%     'timeseriesPopAN', 'timeseriesPopMi');

panelNames = {'B' 'C'};
panelX = {'target onset', 'saccade onset'};
plotPos = {[0.40 0.58 0.55 0.40],[0.40 0.07 0.55 0.40]};
exemplars = [2 22; 2 25; 5 25; 6 14];

%Create a time vector
tVec = cell(1, length(analysisParams.anaWins));
stepSize = analysisParams.binWidth/2.5;
for win_ix = 1:length(analysisParams.anaWins)
    ana_win = analysisParams.anaWins(win_ix);
    edgeStarts = ana_win.winEdges(1):stepSize:ana_win.winEdges(end)-1;
    edgeStops = edgeStarts + analysisParams.binWidth - 1;
    thisEdges = [edgeStarts' edgeStops'];
    while thisEdges(end, 2) > ana_win.winEdges(end) || thisEdges(end, 1) > ana_win.winEdges(end)
        thisEdges = thisEdges(1:end-1, :);
    end
    tVec{win_ix} = mean([thisEdges(:, 1)-1 thisEdges(:, 2)], 2);
end
clear win_ix edgeStarts edgeStops thisEdges

%% Get the statistic and compare it to baseline statistic

if strcmpi(useTsStat, 'F')
    tsStat = timeseriesUnitANP;
    thisBase = cat(1, tuningStatAN{:});
    thisBase = thisBase(:, strcmpi(anaNames, 'baseline'), 1);
elseif strcmpi(useTsStat, 'MI')
    tsStat = timeseriesUnitMiP;
    thisBase = cat(1, tuningStatMI{:});
    thisBase = thisBase(:, strcmpi(anaNames, 'baseline'), 1);
    thisBase = sqrt(thisBase);
end
thisBase = thisBase(keepBool);

pVal = cell(1, length(analysisParams.anaWins));
for win_ix = 1:length(analysisParams.anaWins)
    this_stat = cat(1, tsStat{:, win_ix});
    this_stat = this_stat(keepBool, :, 1);
    if strcmpi(useTsStat, 'MI')
        this_stat = sqrt(this_stat);
    end
    nSteps = size(this_stat, 2);
    pVal{win_ix} = nan(nSteps, 1);
    for step_ix = 1:nSteps
        [~, pVal{win_ix}(step_ix)] = ttest(this_stat(:, step_ix), thisBase,...
            'tail', 'right');
    end
    isSig = pVal{win_ix} < 0.001/nSteps;  %Bonferroni
    sigStart = find(isSig(find(~isSig, 1, 'first'):end), 1, 'first');
    sigStop = find(isSig, 1, 'last');
    
    %Print statement
    if strcmpi(analysisParams.anaWins(win_ix).name, 'timeSeriesTarget')
        lockEventStr = 'target onset';
    elseif strcmpi(analysisParams.anaWins(win_ix).name, 'timeSeriesSaccade')
        lockEventStr = 'saccade onset';
    end
    if win_ix == 2
        fprintf(', and ');
    end
    fprintf('MI was significantly different to baseline from %i to %i msec relative to %s',...
        tVec{win_ix}(sigStart), tVec{win_ix}(sigStop), lockEventStr);
end
fprintf(' (all Bonferroni-corrected p-values < 0.001).\n');

%% Plot the timeseries statistics

% Remember the session id for each neuron.
% This is used for identifying exemplars to plot.
sessId = cell(size(tsStat, 1), 1);
for sess_ix = 1:size(tsStat, 1)
    sessId{sess_ix} = sess_ix*ones(size(tsStat{sess_ix, 1}, 1), 1);
end
sessId = cat(1, sessId{:});
sessId = sessId(keepBool);
clear sess_ix

for win_ix = 1:length(analysisParams.anaWins)
    ana_win = analysisParams.anaWins(win_ix);
    
    this_stat = cat(1, tsStat{:, win_ix});
    this_sig = this_stat(keepBool, :, 2) < 0.001;
    this_stat = this_stat(keepBool, :, 1);  %size nUnits, nTimeSteps
    
    %there are some discontinuities. Stop at the first nan.
    for n_ix = 1:size(this_stat, 1)
        this_stat(n_ix, find(isnan(this_stat(n_ix, :)), 1, 'first'):end) = nan;
    end
    
    if doDeltaMI
        xcut = 0;
        if win_ix == 2
            xcut = -2000;
        end
        this_stat = bsxfun(@minus, this_stat,...
            mean(this_stat(tVec{win_ix}<xcut, :)));
    elseif doRelMI
        this_stat = bsxfun(@rdivide, this_stat, baseMI); 
    end
    
    subplot('Position', plotPos{win_ix});
    hold on
    %scatter plot significant points
    for unit_ix = 1:size(this_stat, 1)
        h = scatter(tVec{win_ix}(this_sig(unit_ix, :)),...
            this_stat(unit_ix, this_sig(unit_ix, :)),...
            24, [0.5 0.5 0.5], 'filled');
        
    end
    plot(tVec{win_ix}, this_stat, 'Color', [0.5 0.5 0.5])
    xlim(ana_win.winEdges);
    if doRelMI
        ylim([0 10]);
    elseif doDeltaMI
        ylim([-0.5 1.5]);
    else
        ylim([0 1.8]);
    end
    xlabel(['Time relative to ' panelX{win_ix} ' (ms)'])
    box off
    set(gca, 'FontSize', fsize_ax)
    set(gca, 'LineWidth', 2)
    set(gca, 'Color', 'none');

    
    % Plot vertical line at time=0
    plot([0 0], get(gca, 'YLim'), 'k--')
    ystr = '';
    if doRelMI
        ystr = 'Relative ';
    end
    if doDeltaMI
        ystr = [ystr '\Delta '];
    end
    ystr = [ystr 'Mutual Information'];
    ylabel(ystr);
    
    xl = get(gca, 'XLim');
    yl = get(gca, 'YLim');
    h = text(min(xl) - diff(xl)/20, max(yl), panelNames{win_ix});
    set(h, 'FontSize', fsize_lab);
    
    % Plot exemplar
    h = [];
    leg_names = cell(1, size(exemplars, 1));
    for ex_ix = 1:size(exemplars, 1)
        unit_ix = find(sessId == exemplars(ex_ix, 1), exemplars(ex_ix, 2), 'first');
        unit_ix = unit_ix(end);
        
        x_h = plot(tVec{win_ix}, this_stat(unit_ix, :)', 'LineWidth', 3);
        h = [h, x_h];
        leg_names{ex_ix} = [num2str(exemplars(ex_ix, 1)) '.' num2str(exemplars(ex_ix,2))];
        %Scatter-plot significant points
        s_h = scatter(tVec{win_ix}(this_sig(unit_ix, :)),...
            this_stat(unit_ix, this_sig(unit_ix, :))',...
            100, x_h.Color, 'filled');
        s_h.MarkerEdgeColor = 'k';
        
    end
    if win_ix == 1
        legend(h, leg_names);
        legend boxoff
    end
    
    % Inidicate significant time-steps
    x = tVec{win_ix};
    y = (max(yl)-1.1*diff(yl)/10) * ones(size(x));
    y(pVal{win_ix} > 0.001/nSteps) = nan;
    plot(x, y, 'k', 'LineWidth', 2);
    
    % Indicate analysis epochs.
    if win_ix == 1
        epoch_names = {'baseline' 'cue', 'delay'};
        per_cols = {'w' 'g' 'b'};
    elseif win_ix == 2
        epoch_names = {'resp'};
        per_cols = {'r'};
    end
    for ep_ix = 1:length(epoch_names)
        per_bool = strcmpi({analysisEpochs.name}, epoch_names{ep_ix});
        x = repmat(analysisEpochs(per_bool).winEdges, 2, 1);
        y = [max(yl) max(yl)-diff(yl)/10 max(yl)-diff(yl)/10 max(yl)];
        hp = patch(x(:), y(:), per_cols{ep_ix});%, 'EdgeColor', 'none');
        if per_cols{ep_ix} == 'b'
            set(hp, 'FaceAlpha', 0.5);
        end
        ht = text(mean(analysisEpochs(per_bool).winEdges), mean(y(:)), epoch_names{ep_ix}, 'FontSize', fsize_ax);
        set(ht, 'HorizontalAlignment', 'center');
    end
    hold off
end

%%
set(miFig, 'Color', 'none');
savefig(miFig, fullfile(paths.results, 'Figures', 'mi_timeseries'));