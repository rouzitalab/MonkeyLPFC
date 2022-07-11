clear
clc
path = 'G:\Projects\gpfa_v0203\data';
data_list = dir(fullfile(path, '*.mat'));
data_path = append(path, '\',data_list(1).name);
data1 = load(data_path);

lstm_out = abs(data1.lstm);
label = data1.label;
y = data1.y;
spikes = data1.spikes;
trial = data1.trial;

runIdx = 40;
for i = 1 : size(lstm_out,1)
    dat(i).trialId = trial(i);
    dat(i).spikes = squeeze(spikes(i,:,:));
end

% Select method to extract neural trajectories:
% 'gpfa' -- Gaussian-process factor analysis
% 'fa'   -- Smooth and factor analysis
% 'ppca' -- Smooth and probabilistic principal components analysis
% 'pca'  -- Smooth and principal components analysis
method = 'gpfa';

% Select number of latent dimensions
xDim = 8;
% NOTE: The optimal dimensionality should be found using 
%       cross-validation (Section 2) below.

% If using a two-stage method ('fa', 'ppca', or 'pca'), select
% standard deviation (in msec) of Gaussian smoothing kernel.
kernSD = 60;
% NOTE: The optimal kernel width should be found using 
%       cross-validation (Section 2) below.

% Extract neural trajectories
result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,... 
                    'binWidth', kernSD);
% NOTE: This function does most of the heavy lifting.

% Orthonormalize neural trajectories
[estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
% NOTE: The importance of orthnormalization is described on 
%       pp.621-622 of Yu et al., J Neurophysiol, 2009.
%% Smoothing with original points of interest
% for i = 1 : length(seqTrain)
%     tmp = seqTrain(i).xorth;
%     tmp2 = smoother(tmp,1000, 100, 'causal', true);
%     tmp = seqTrain(i).xorth;
%     seqTrain(i).xorth = tmp2;
%     seqTrain(i).xorth(:,6) = tmp(:,6);
%     seqTrain(i).xorth(:,9) = tmp(:,9);
%     seqTrain(i).xorth(:,19) = tmp(:,19);
% end
%% Averaging points of interest among the same targets
% u = idx(find(direc(idx)==up));
% d = idx(find(direc(idx)==down));
% cold = 0; colu = 0; gou = 0; god = 0;
% for i = 1:length(seqTrain)
% if ismember(seqTrain(i).trialId, u), colu = colu + seqTrain(i).xorth(:,9); gou = gou + seqTrain(i).xorth(:,19);
% else, cold = cold + seqTrain(i).xorth(:,9); god = god + seqTrain(i).xorth(:,19); end
% end
% colu = colu/length(u);gou = gou/length(u);
% cold = cold/length(d);god = god/length(d);
% for i = 1:length(seqTrain)
% if ismember(seqTrain(i).trialId, u), seqTrain(i).xorth(:,9)=colu; seqTrain(i).xorth(:,19)=gou;
% else, seqTrain(i).xorth(:,9)=cold; seqTrain(i).xorth(:,19)=god; end
% end
%% Averaging points of interest among the same rule
% cr = 0; gor = 0;
% for i = 1:length(seqTrain)
%     cr = cr + seqTrain(i).xorth(:,9); gor = gor + seqTrain(i).xorth(:,19);
% end
% cr = cr / length(seqTrain); gor = gor / length(seqTrain);
% for i = 1:length(seqTrain)
%     seqTrain(i).xorth(:,9)=cr; seqTrain(i).xorth(:,19)=gor;
% end
%% Plot neural trajectories in 2D and 3D spaces
plot2D(seqTrain, label,'xorth', 'dimsToPlot', 1:2);
plot3D(seqTrain, label,'xorth', 'dimsToPlot', 1:3);