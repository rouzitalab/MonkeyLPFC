% =====
% as_gpfa_example3
% example 3 does analysis on actual array data generated from JerryLee etc
% tt is for t-test as opposed to decoder
% TIPS
% =====
% - For exploratory analysis using GPFA, we often run only Section 1
%   below, and not Section 2 (which finds the optimal laten
%   dimensionality).  This can provide a substantial savings in running
%   time, since running Section 2 takes roughly K times as long as
%   Section 1, where K is the number of cross-validation folds.  As long
%   as we use a latent dimensionality that is 'large enough' in Section 1,
%   we can roughly estimate the latent dimensionality by looking at
%   the plot produced by plotEachDimVsTime.m.  The optimal latent
%   dimensionality is approximately the number of top dimensions that
%   have 'meaningful' temporal structure.  For visualization purposes,
%   this rough dimensionality estimate is usually sufficient.
%
% - For exploratory analysis with the two-stage methods, we MUST run
%   Section 2 to obtain the optimal smoothing kernel width.  There is
%   no easy way estimate the optimal smoothing kernel width from the
%   results of Section 1.

% ===========================================
% 1) Basic extraction of neural trajectories
% ===========================================
% load('mat_sample/sample_dat');

% Results will be saved in mat_results/runXXX/, where XXX is runIdx.
% Use a new runIdx for each dataset.

% Go to em.m to change number of EM iterations to 500
% model.Mu   : This variable gives the cluster-dependent regression
%   coefficients. .Mu(:,k,d) gives the coefficients for the k-th cluster
%   in the d-th dimension; note that .Mu(1,k,d) gives the y-intercept.
## first do this in 2 (make 2tt), then go through each trajectory individually (ie with for loop) and fit 
whatever to it, then make a big matrix with fitted parameters for stat testing

runIdx = 979;
analysisBlock=4; %[1 8];%4; these are groups of analysis Periods
analysisGroup=[];
%analysisGroup{1}=[1 2 3 4 5 6 7 8];  %these are in analysis Periods
%analysisGroup{2}=[57 58 59 60 61 62 63 64];
% analysisGroup{1}=[ 3  6]+1;  %these are in analysis Periods
% analysisGroup{2}=[ 59  62 ]+1;
excludeNoiseUnits = 1;
ops.K=8;  %# of allowable clusters
%filename =  'sra3_2_j_038_00+02.resmat' %'sra3_2_j_039_00+01.resmat'%'
filename =  'sra3_2_j_039_00+01.resmat'
%apFileName = 'anaParams20\';
%apFileName = 'anaP10_saccadehit';
respath = '/Users/adam/Desktop/Matlab add ons/resmat repository/';  % for macintosh
%apFileName = 'anaParams20/';
apFileName = 'anaP10_saccadehit/';
%respath = 'c:\Documents and Settings\Administrator\Desktop\resmat repository\';

path = strcat(respath,apFileName);
cd (strcat(path))

% next line changed from as_convertres2dat1..
%dat1=as_convertres2dat2(filename,analysisBlock,excludeNoiseUnits);
dat1=as_convertres2dat2(filename,analysisBlock,excludeNoiseUnits,analysisGroup);
% dat1 = as_SimSaccadeExp(80,0.8);
% step one, order the trials
ndir=8;
index1=zeros(size(dat1,2),2);




% Select method to extract neural trajectories:
% 'gpfa' -- Gaussian-process factor analysis
% 'fa'   -- Smooth and factor analysis
% 'ppca' -- Smooth and probabilistic principal components analysis
% 'pca'  -- Smooth and principal components analysis
method = 'gpfa'; %'gpfa';
if ~ exist('analysisGroup') || isempty(analysisGroup)
    for i = 1:size(dat1,2)
        index1(i,:)=[dat1(i).dir dat1(i).trialId];
    end
    [Y,I]=sort(index1(:,1));
    
    
    for j=1:size(dat1,2)
        dat(j)=dat1(I(j));
        dat(j).trialId = j;
        
    end
    
    origidx = Y;
    for i =1:ndir
        origGroup{i}=i*ones(sum(Y==i),1);
    end
    for iLookup = 1: size(dat,2)
        trialDirLookup(iLookup,:)=[dat(iLookup).trialId dat(iLookup).dir];
    end
else
    for i = 1:size(dat1,2)
        index1(i,:)=[dat1(i).analysisPeriod dat1(i).trialId];
    end
    [Y,I]=sort(index1(:,1));
    
    
    for j=1:size(dat1,2)
        dat(j)=dat1(I(j));
        dat(j).trialId = j;
        
    end
    
    origidx = Y;
    for i =1:ndir
        origGroup{i}=i*ones(sum(Y==i),1);
    end
    for iLookup = 1: size(dat,2)
        trialDirLookup(iLookup,:)=[dat(iLookup).trialId dat(iLookup).analysisPeriod];
    end
    
end
% Select number of latent dimensions
xDim = 8; %coincidence ajs
% NOTE: The optimal dimensionality should be found using
%       cross-validation (Section 2) below.

% If using a two-stage method ('fa', 'ppca', or 'pca'), select
% standard deviation (in msec) of Gaussian smoothing kernel.
kernSD = 60; %30;
% NOTE: The optimal kernel width should be found using
%       cross-validation (Section 2) below.

% Extract neural trajectories
result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,...
    'kernSDList', kernSD);
% NOTE: This function does most of the heavy lifting.
%
% Orthonormalize neural trajectories
[estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
% NOTE: The importance of orthnormalization is described on
%       pp.621-622 of Yu et al., J Neurophysiol, 2009.

% Plot neural trajectories in 3D space
as_plot3D(trialDirLookup, seqTrain, 'xorth', 'dimsToPlot', 1:3);
as_plot2D(trialDirLookup, seqTrain, 'xorth', 'dimsToPlot', [1 2]); %1:2

%for i=1:size(dat,2)
%   Seq2DatIdx(i)=seqTrain(i).trialId;
%end

%for i = 1:size(dat,2)
%   Trajs.Y{i}=[seqTrain(Seq2DatIdx(i)).xorth(1,:);seqTrain(Seq2DatIdx(i)).xorth(2,:)]';
%end
%
for i=1:size(seqTrain,2)
    seqTrain1(seqTrain(i).trialId)=seqTrain(i);
end

%%
%   for i = 1:size(dat,2)
%       Trajs.Y{i}=[seqTrain1(i).xorth(3,:);seqTrain1(i).xorth(2,:)]';
%   end
 for i = 1:size(dat,2)
     Trajs.Y{i}=[seqTrain1(i).xorth(1,:);seqTrain1(i).xorth(2,:);seqTrain1(i).xorth(3,:)]'; % TRY FITTING TO MULTIDIMENSION
 
 end
%%
ops.K=ops.K;

ops.zero='none';
ops.order=1;

ops.method = 'srm'; %'lrm' 'gmix' 'srm' ListModels()
ops.IterLimit=200;
ops.NumEMStarts=8;
ops.ShowGraphics	=0;
model = curve_clust(Trajs,ops); % THIS IS THE CURVE CLUSTERING WORK HORSE
showmodel(model,Trajs);
AIC=2*model.K - 2*model.TrainLhood_ppt;
idx=model.C;


%%
% x = meas;
for i = 1:ops.K
    unorderedGroupEst{i} = origidx(idx == i)';
end

% [orderedGroupEst,map, cm] = as_classify3(origGroup,unorderedGroupEst,1);
[orderedGroupEst,map, cm] = as_classify2(origidx,unorderedGroupEst);

H = as_infoClust1(cm)
Pc = sum(diag(cm)) / sum(sum(cm))
 
 
 
 
 
 
 
 %plot(seqTrain, 'xorth', 'dimsToPlot', 1:2);
% NOTES:
% - This figure shows the time-evolution of neural population
%   activity on a single-trial basis.  Each trajectory is extracted from
%   the activity of all units on a single trial.
% - This particular example is based on multi-electrode recordings
%   in premotor and motor cortices within a 400 ms period starting 300 ms 
%   before movement onset.  The extracted trajectories appear to
%   follow the same general path, but there are clear trial-to-trial
%   differences that can be related to the physical arm movement. 
% - Analogous to Figure 8 in Yu et al., J Neurophysiol, 2009.
% WARNING:
% - If the optimal dimensionality (as assessed by cross-validation in 
%   Section 2) is greater than 3, then this plot may mask important 
%   features of the neural trajectories in the dimensions not plotted.  
%   This motivates looking at the next plot, which shows all latent 
%   dimensions.

% Plot each dimension of neural trajectories versus time
% plotEachDimVsTime(seqTrain, 'xorth', result.binWidth);
% NOTES:
% - These are the same neural trajectories as in the previous figure.
%   The advantage of this figure is that we can see all latent
%   dimensions (one per panel), not just three selected dimensions.  
%   As with the previous figure, each trajectory is extracted from the 
%   population activity on a single trial.  The activity of each unit 
%   is some linear combination of each of the panels.  The panels are
%   ordered, starting with the dimension of greatest covariance
%   (in the case of 'gpfa' and 'fa') or variance (in the case of
%   'ppca' and 'pca').
% - From this figure, we can roughly estimate the optimal
%   dimensionality by counting the number of top dimensions that have
%   'meaningful' temporal structure.   In this example, the optimal 
%   dimensionality appears to be about 5.  This can be assessed
%   quantitatively using cross-validation in Section 2.
% - Analogous to Figure 7 in Yu et al., J Neurophysiol, 2009.

fprintf('\n');
fprintf('Basic extraction and plotting of neural trajectories is complete.\n');
fprintf('Press any key to start cross-validation...\n');
fprintf('[Depending on the dataset, this can take many minutes to hours.]\n');
% pause;


% ========================================================
% 2) Full cross-validation to find:
%  - optimal state dimensionality for all methods
%  - optimal smoothing kernel width for two-stage methods
% ========================================================

if 0
% Select number of cross-validation folds
numFolds = 4;

% Perform cross-validation for different state dimensionalities.
% Results are saved in mat_results/runXXX/, where XXX is runIdx.
for xDim = [2 5 8]
  neuralTraj(runIdx, dat, 'method',  'pca', 'xDim', xDim, 'numFolds', numFolds);
  neuralTraj(runIdx, dat, 'method', 'ppca', 'xDim', xDim, 'numFolds', numFolds);
  neuralTraj(runIdx, dat, 'method',   'fa', 'xDim', xDim, 'numFolds', numFolds);
  neuralTraj(runIdx, dat, 'method', 'gpfa', 'xDim', xDim, 'numFolds', numFolds);
end
fprintf('\n');
% NOTES:
% - These function calls are computationally demanding.  Cross-validation 
%   takes a long time because a separate model has to be fit for each 
%   state dimensionality and each cross-validation fold.

% Plot prediction error versus state dimensionality.
% Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
kernSD = 30; % select kernSD for two-stage methods
plotPredErrorVsDim(runIdx, kernSD);
% NOTES:
% - Using this figure, we i) compare the performance (i.e,,
%   predictive ability) of different methods for extracting neural
%   trajectories, and ii) find the optimal latent dimensionality for
%   each method.  The optimal dimensionality is that which gives the
%   lowest prediction error.  For the two-stage methods, the latent
%   dimensionality and smoothing kernel width must be jointly
%   optimized, which requires looking at the next figure.
% - In this particular example, the optimal dimensionality is 5. This
%   implies that, even though the raw data are evolving in a
%   53-dimensional space (i.e., there are 53 units), the system
%   appears to be using only 5 degrees of freedom due to firing rate
%   correlations across the neural population.
% - Analogous to Figure 5A in Yu et al., J Neurophysiol, 2009.

% Plot prediction error versus kernelSD.
% Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
xDim = 5; % select state dimensionality
plotPredErrorVsKernSD(runIdx, xDim);
% NOTES:
% - This figure is used to find the optimal smoothing kernel for the
%   two-stage methods.  The same smoothing kernel is used for all units.
% - In this particular example, the optimal standard deviation of a
%   Gaussian smoothing kernel with FA is 30 ms.
% - Analogous to Figures 5B and 5C in Yu et al., J Neurophysiol, 2009.
end
