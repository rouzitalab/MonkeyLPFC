function reducedD = reduceDimsAndTransformWithC(D, C, dhParams, dhRedParams, mode)

%Convert testing dat to transformable format
%using reducedims with dhRedMode == -1
[newD, ~, ~, ~] = reducedims(D, -1, size(D(1).data, 1), dhParams);
[newD.y] = deal(newD.data);  % Copy data into y.
t_length = nan(size(newD));
for tr_ix = 1:numel(newD)
    t_length(tr_ix) = size(newD(tr_ix).data, 2);
end
t_length = num2cell(t_length);

tId = num2cell(trialIds);
[newD.trialId] = tId{:};
[newD.T] = t_length{:};
clear tr_ix t_length

if any(strcmpi({'pca', 'fa', 'lda'}, mode))
    offsets = mean([newD.data], 2);
    %offsets = zeros(size([testD.data], 1), 1);
    for tr_ix = 1:length(newD)
        newD(tr_ix).data = C' * bsxfun(@minus, newD(tr_ix).data, offsets);
    end; clear tr_ix offsets
    reducedD = newD;
elseif strcmpi('gpfa', mode)
    reducedD = exactInferenceWithLL(newD, dhRedParams);
    [Xorth, Corth] = orthogonalize([reducedD.xsm], dhRedParams.C);
    reducedD = segmentByTrial(reducedD, Xorth, 'data');
    reducedD = rmfield(reducedD, {'Vsm', 'VsmGP', 'xsm'});
elseif any(strcmpi({'dPCA', 'canon'}, mode))
    reducedD = centerAndTransformStructArray(newD, C);
end