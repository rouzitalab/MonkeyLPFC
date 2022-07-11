function mystruct = centerAndTransformStructArray(mystruct, W)
%mystruct = centerandtransformstructarray(mystruct)
%mystruct is a struct array with each element is a trial containing a
% .data field. The data field is a matrix of shape m_sources x n_times

[mSources, ~] = size(mystruct(1).data);
nTrials = length(mystruct);
% Get the mean across trials
acrossTrialSum = zeros(mSources, 1);
sampCount = nan(1, nTrials);
for tr_ix = 1:nTrials
    acrossTrialSum = acrossTrialSum + sum(mystruct(tr_ix).data, 2);
    sampCount(tr_ix) = size(mystruct(tr_ix).data, 2);
end
acrossTrialMeans = acrossTrialSum ./ sum(sampCount);
for tr_ix = 1:nTrials
    tmp = bsxfun(@minus, mystruct(tr_ix).data, acrossTrialMeans);
    mystruct(tr_ix).data = W' * tmp;
end