function D = cDat2DataHigh(trials)

D = struct('data', [], 'epochStarts', [], 'epochColors', [], 'condition', {trials.newClassStr}, 'type', 'traj');

% We need to map the class to a color. For now, assume ordering newClass
% will work.
uqClass = unique([trials.newClass]);
nClasses = length(uqClass);
[~, classIx] = ismember([trials.newClass], uqClass);
cmap = parula(2*nClasses);

for tr_ix = 1:length(D)
    tr = trials(tr_ix);
    D(tr_ix).data = full(tr.raster(tr.timeBool, :)');
    trT = tr.tVec(tr.timeBool);
    epoch_starts = find(abs(trT) == min(abs(trT)), 1, 'first');
    cmap_offs = -1;
    if ~isempty(tr.avoidTime) && ~isnan(tr.avoidTime)...
            && round(tr.avoidTime - tr.eventTime) + epoch_starts(1) < sum(tr.timeBool)
        epoch_starts = [epoch_starts round(tr.avoidTime - tr.eventTime) + epoch_starts(1)];
        cmap_offs = -1:0;
    end
    this_cmap = cmap(2*classIx(tr_ix) + cmap_offs, :);
    D(tr_ix).epochStarts = epoch_starts;
    D(tr_ix).epochColors = this_cmap;
end