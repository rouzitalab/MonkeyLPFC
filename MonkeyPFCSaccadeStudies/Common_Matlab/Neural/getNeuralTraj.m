function [estParams, seqTrain, seqTest] = getNeuralTraj(seq, testMask, varargin )
%[estParams, seqTrain, seqTest] = getNeuralTraj(seq, testMask, 'param', value )
%testMask is a boolean vector equal to length seq
%validParamNames = {'winEdges' 'binWidth' 'runIdx' 'method' 'xDim' 'kernSD'}

%% Parameters
%Default params.
params.winEdges = [-300 100];
params.binWidth = 20;
params.runIdx = 1;% Results will be saved in mat_results/runXXX/, where XXX is runIdx.
params.method = 'gpfa'; %1-stage: 'gpfa'; 2-stage: 'fa' 'ppca' 'pca'
params.xDim = 8;% Select number of latent dimensions, optimal value found using cross-validation.
params.kernSD = 30; %msec smoothing for 2-stage, optimal value found using cross-validation.
%Consolidate with varargin
params = varg2params(varargin, params, {'winEdges' 'binWidth' 'runIdx' 'method' 'xDim' 'kernSD'});

%Will be saved into result.
winEdges = params.winEdges;
binWidth = params.binWidth;
runIdx = params.runIdx;
method = params.method;
xDim = params.xDim;
kernSD = params.kernSD;
cvf = 0;
extraOpts = {};


%% Housekeeping

if ~isdir(fullfile(pwd, 'mat_results'))
    mkdir(fullfile(pwd, 'mat_results'));
end

runDir = fullfile(pwd,'mat_results',['run' sprintf('%03d',params.runIdx)]);

if isdir(runDir)
    rmdir(runDir,'s');
end
mkdir(runDir);

fname = sprintf('%s/%s_xDim%02d', runDir, params.method, params.xDim);

if isempty(seq)
    fprintf('Error: No valid trials.  Exiting.\n');
    result = [];
    return;
end

%% Extract neural trajectories

seqTrain = seq(~testMask);
seqTest = seq(testMask);

% Remove inactive units based on training set
hasSpikesBool = (mean([seqTrain.y], 2) ~= 0);

for n = 1:length(seqTrain)
  seqTrain(n).y = seqTrain(n).y(hasSpikesBool,:);
end
for n = 1:length(seqTest)
  seqTest(n).y = seqTest(n).y(hasSpikesBool,:);
end

% Check if training data covariance is full rank
yAll = [seqTrain.y];
yDim  = size(yAll, 1);
    
if rank(cov(yAll')) < yDim
  fprintf('ERROR: Observation covariance matrix is rank deficient.\n');
  fprintf('Possible causes: repeated units, not enough observations.\n');
  fprintf('Exiting...\n');
  return
end

fprintf('Number of training trials: %d\n', length(seqTrain));
fprintf('Number of test trials: %d\n', length(seqTest));
fprintf('Latent space dimensionality: %d\n', params.xDim);
fprintf('Observation dimensionality: %d\n', sum(hasSpikesBool));

% The following does the heavy lifting.
if isequal(params.method, 'gpfa')
  gpfaEngine(seqTrain, seqTest, fname,... 
  'xDim', params.xDim, 'binWidth', params.binWidth); %startTau, startEps??

elseif ismember(params.method, {'fa', 'ppca', 'pca'})
  twoStageEngine(seqTrain, seqTest, fname,... 
  'typ', params.method, 'xDim', params.xDim, 'binWidth', params.binWidth);  %kernSDList
end

if exist([fname '.mat'], 'file')
  save(fname, 'method', 'cvf', 'hasSpikesBool', 'extraOpts', '-append');
end

result = [];  
if exist([fname '.mat'], 'file')
    result = load(fname);
end

%% Orthonormalize neural trajectories
[estParams, seqTrain, seqTest] = postprocess(result, 'kernSD', params.kernSD);
% NOTE: The importance of orthnormalization is described on 
%       pp.621-622 of Yu et al., J Neurophysiol, 2009.

% 
% %% Plot neural trajectories in 3D space
% plot3D(seqTrain, 'xorth', 'dimsToPlot', 1:3);
% % NOTES:
% % - This figure shows the time-evolution of neural population
% %   activity on a single-trial basis.  Each trajectory is extracted from
% %   the activity of all units on a single trial.
% % - This particular example is based on multi-electrode recordings
% %   in premotor and motor cortices within a 400 ms period starting 300 ms 
% %   before movement onset.  The extracted trajectories appear to
% %   follow the same general path, but there are clear trial-to-trial
% %   differences that can be related to the physical arm movement. 
% % - Analogous to Figure 8 in Yu et al., J Neurophysiol, 2009.
% % WARNING:
% % - If the optimal dimensionality (as assessed by cross-validation in 
% %   Section 2) is greater than 3, then this plot may mask important 
% %   features of the neural trajectories in the dimensions not plotted.  
% %   This motivates looking at the next plot, which shows all latent 
% %   dimensions.
% 
% % Plot each dimension of neural trajectories versus time
% plotEachDimVsTime(seqTrain, 'xorth', result.binWidth);
% % NOTES:
% % - These are the same neural trajectories as in the previous figure.
% %   The advantage of this figure is that we can see all latent
% %   dimensions (one per panel), not just three selected dimensions.  
% %   As with the previous figure, each trajectory is extracted from the 
% %   population activity on a single trial.  The activity of each unit 
% %   is some linear combination of each of the panels.  The panels are
% %   ordered, starting with the dimension of greatest covariance
% %   (in the case of 'gpfa' and 'fa') or variance (in the case of
% %   'ppca' and 'pca').
% % - From this figure, we can roughly estimate the optimal
% %   dimensionality by counting the number of top dimensions that have
% %   'meaningful' temporal structure.   In this example, the optimal 
% %   dimensionality appears to be about 5.  This can be assessed
% %   quantitatively using cross-validation in Section 2.
% % - Analogous to Figure 7 in Yu et al., J Neurophysiol, 2009.
% 
% 
% % ========================================================
% % 2) Full cross-validation to find:
% %  - optimal state dimensionality for all methods
% %  - optimal smoothing kernel width for two-stage methods
% % ========================================================
% 
% % Select number of cross-validation folds
% numFolds = 4;
% 
% % Perform cross-validation for different state dimensionalities.
% % Results are saved in mat_results/runXXX/, where XXX is runIdx.
% for xDim = [2 5 8 15]
%   neuralTraj(params.runIdx, dat, 'method',  'pca', 'xDim', xDim, 'numFolds', numFolds);
%   neuralTraj(params.runIdx, dat, 'method', 'ppca', 'xDim', xDim, 'numFolds', numFolds);
%   neuralTraj(params.runIdx, dat, 'method',   'fa', 'xDim', xDim, 'numFolds', numFolds);
%   neuralTraj(params.runIdx, dat, 'method', 'gpfa', 'xDim', xDim, 'numFolds', numFolds);
% end
% fprintf('\n');
% % NOTES:
% % - These function calls are computationally demanding.  Cross-validation 
% %   takes a long time because a separate model has to be fit for each 
% %   state dimensionality and each cross-validation fold.
% 
% % Plot prediction error versus state dimensionality.
% % Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
% kernSD = 30; % select kernSD for two-stage methods
% plotPredErrorVsDim(params.runIdx, kernSD);
% % NOTES:
% % - Using this figure, we i) compare the performance (i.e,,
% %   predictive ability) of different methods for extracting neural
% %   trajectories, and ii) find the optimal latent dimensionality for
% %   each method.  The optimal dimensionality is that which gives the
% %   lowest prediction error.  For the two-stage methods, the latent
% %   dimensionality and smoothing kernel width must be jointly
% %   optimized, which requires looking at the next figure.
% % - In this particular example, the optimal dimensionality is 5. This
% %   implies that, even though the raw data are evolving in a
% %   53-dimensional space (i.e., there are 53 units), the system
% %   appears to be using only 5 degrees of freedom due to firing rate
% %   correlations across the neural population.
% % - Analogous to Figure 5A in Yu et al., J Neurophysiol, 2009.
% 
% % Plot prediction error versus kernelSD.
% % Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
% xDim = 5; % select state dimensionality
% plotPredErrorVsKernSD(params.runIdx, xDim);
% % NOTES:
% % - This figure is used to find the optimal smoothing kernel for the
% %   two-stage methods.  The same smoothing kernel is used for all units.
% % - In this particular example, the optimal standard deviation of a
% %   Gaussian smoothing kernel with FA is 30 ms.
% % - Analogous to Figures 5B and 5C in Yu et al., J Neurophysiol, 2009.

end

