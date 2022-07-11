function [DataHigh, W, V, whichMarg, params] = dPCAofDataHigh(DataHigh, varargin)
%[DataHigh, W, V, whichMarg, params] = dPCAofDataHigh(DataHigh, varargin)

%% Parameters
%Default params.
params.optimLambda = 0.0025;
params.doPlot = false;
params.margCombs = {{1, [1 2]}, {2}};
params.margNames = {'Location-Time', 'Time'};

%Consolidate with varargin
params = varg2params(varargin, params, {'optimLambda', 'doPlot',...
    'margCombs', 'margNames', 'maxDims'});

%% dPCA expects a 3D matrix of size Neurons x Class x Time
[nNeurons, nTimes] = size(DataHigh(1).data);
[uqClass, ~, classIx] = unique({DataHigh.condition});
classCount = hist(classIx, unique(classIx));
maxTrials = max(classCount);

dpcaX = nan(nNeurons, length(uqClass), nTimes, maxTrials);
for cl_ix = 1:length(uqClass)
    trBool = strcmpi({DataHigh.condition}, uqClass{cl_ix});
    dpcaX(:, cl_ix, :, 1:sum(trBool)) = cat(3, DataHigh(trBool).data);
end
dpcaXAvg = nanmean(dpcaX, 4);

%% Run dPCA
trialNum = ones(nNeurons, 1) * classCount;
if params.optimLambda < 0
    params.optimLambda = dpca_optimizeLambda(dpcaXAvg, dpcaX, trialNum, ...
        'combinedParams', params.margCombs, ...
        'numComps', min([nNeurons nTimes 15]),...
        'numRep', 100, ...  % increase this number to ~10 for better accuracy
        'filename', 'tmp_optimalLambdas.mat');
end

[W,V,whichMarg] = dpca(dpcaXAvg, min([nNeurons nTimes 15]), ...
    'combinedParams', params.margCombs, ...
    'lambda', params.optimLambda);

if params.doPlot
    explVar = dpca_explainedVariance(dpcaXAvg, W, V, ...
        'combinedParams', params.margCombs);
    dpca_plot(dpcaXAvg, W, V, @dpca_plot_default, ...
        'explainedVar', explVar, ...
        'marginalizationNames', params.margNames, ...
        'whichMarg', whichMarg,                 ...
        'timeMarginalization', 2, ...
        'legendSubplot', 12);
end