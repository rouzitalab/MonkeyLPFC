 %% Spikes time
% for m74: 200 -> target onset
%           450 -> color cue
%           1450 -> delay
%           1700 -> response
% for j51: 600 -> target onset
%           850 -> color cue
%           1850 -> delay
%           2100 -> response

%% Loading the data

for sess = 1 : 7
    for r = 1:2
        for direc = 0:7
            clearvars -except 'sess' 'r' 'direc'
            %% Load data
            dataset = ["sra3_1_j_050_00_v1_segmented.mat","sra3_1_j_051_00_v1_segmented.mat","sra3_1_j_052_00_v1_segmented.mat","sra3_1_m_077_0001_v1_segmented.mat","sra3_1_m_081_0001_v1_segmented.mat","sra3_1_m_082_0001_v1_segmented.mat","sra3_1_m_083_0001_v1_segmented.mat"];
            path = 'G:\Projects\gpfa_v0203\data';
            data_list = dir(fullfile(path, '*.mat'));
            runid = 108 + 16*(sess-1) + 8*(r-1) + direc + 1;
            data_path = append(path, '\',dataset(sess));
            data1 = load(data_path);
            % data2 = load(fullfile(data_path),'-mat');

            % outcome = data1.outcome;
            % flag = find(outcome~=-1);
            spike = data1.spikes;
            % spike(isnan(spike))=0;
            spikes = permute(spike,[1,3,2]);
            % spikes = int8(spikes);
            saccades = data1.saccade;
            targets = data1.target;
            colors = cell2mat(data1.color);
            learned = data1.learned + 1;
            unlearned = data1.unlearned + 1;
            %% change every iteration
            runIdx = runid;
            if r == 1
                rule = learned;
            else
                rule = unlearned;
            end
            % dir = 6;
            %% Learned and unlearned trials
            % [learned,unlearned, border]=monkey_performance(data1);
            % ups = [1,1];
            % downs = [4,5,6,7];
            % colors = ['g','r'];
            % for run = 1:2
            % spikes = data1.spikes(flag,:,:);
            % saccade = data1.saccade(flag);
            % target = data1.target(flag);
            % color = cell2mat(data1.color(flag));
            %% Change per run

            % trials = learned;
            % trials = unlearned;
            % up = ups(run);
            % down = downs(run);
            % colr = colors(run);
            % direc = saccade;
            % direc = target;
            % direction = direc(trials);
            % color = color(trials);
            %% Devide trials by targets
            % idx = trials(find((direction==up)|(direction==down)));
            %% Devide trials by rule
            % idx = trials(find((direction==up)&(color==colr)));
            %% Devide trials by block
            % tmp = find(trials>border(1));
            % tmp2 = trials(tmp);
            % tmp = find(b<border(2));
            % idx = tmp2(tmp);
            %% Select spike trials 
            % spikes = spikes(idx,:,:);
            % spikes = spikes(learned,:,:);
            % spikes(isnan(spikes))=0;
            %% Downsampling without loosing spikes
            % k = find(spikes);
            % k = floor(k/5);
            % if k(1)==0, k(1)=1; end
            % spikeTrain = zeros(size(spikes,1),floor(size(spikes,2)/5)+1, size(spikes,3));
            % % spikeTrain = zeros(size(spikes));
            % spikeTrain(k)=1;
            % spikeTrain = spikes;
            %% Data structure for the method

            idx = targets(rule)==direc;
            trial = rule(idx);
            tmp = spikes(rule,:,:);
            spikeTrain = tmp(idx,:,:);
            tmp = colors(rule);
            color = tmp(idx);
            for i = 1 : length(trial)
            %     dat(i).trialId = idx(i);
                dat(i).trialId = trial(i);
                dat(i).spikes = squeeze(spikeTrain(i,:,:));
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
            kernSD = 20;
            % NOTE: The optimal kernel width should be found using 
            %       cross-validation (Section 2) below.

            % Extract neural trajectories
            result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,... 
                                'kernSDList', kernSD);
            % NOTE: This function does most of the heavy lifting.

            % Orthonormalize neural trajectories
            [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
            % NOTE: The importance of orthnormalization is described on 
            %       pp.621-622 of Yu et al., J Neurophysiol, 2009.
            % v = ones(1,5);
            % for i = 1 : length(seqTrain)
            %     strc(i).data = seqTrain(i).xorth;
            % end
            % spike2d = zeros(size(spikeTrain,1),size(spikeTrain,3)*size(spikeTrain,2));
            % fr = zeros(size(spikeTrain,1),size(spikeTrain,3)*size(spikeTrain,2));
            % for i = 1 : size(spikeTrain,1)
            %     for k = 1 : size(spikeTrain,2)
            %         spike2d(i, (k-1)*size(spikeTrain,3)+1:k*size(spikeTrain,3)) = spikeTrain(i,k,:);
            %     end
            % end
            % for i = 1 : size(spike2d,1)
            %     for j = 25 : size(spike2d,2)
            %         fr(i,j)=sum(spike2d(i,j-9:j));
            %     end
            % end
            % fr = fr(learned,:);
            % spike2d = spike2d(learned,:);
            % sac = saccades(learned)';
            % [~,~,stats] = manova1(fr,sac);
            % W = stats.eigenvec(:,1:8);
            % trajectories = centerAndTransformStructArray(strc,W);

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
            %% Plot title
            set = extractBefore(dataset(sess),strlength(dataset(sess))-3);
            if r == 1
                rul = 'learned';
            else
                rul = 'unlearned';
            end
            titl = append('Direction = ', int2str(direc), ' for ', rul, ' trials in ', set);

            %% Plot neural trajectories in 2D and 3D spaces
            plot2D(seqTrain, color,'xorth',titl, 'dimsToPlot', 1:2);
        end
    end
end
% plot3D(seqTrain, saccades,'xorth', 'dimsToPlot', 1:3);
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
% plotEachDimVsTime(seqTrain, 'xorth',saccades, result.binWidth);
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
% end
fprintf('\n');
fprintf('Basic extraction and plotting of neural trajectories is complete.\n');
fprintf('Press any key to start cross-validation...\n');
fprintf('[Depending on the dataset, this can take many minutes to hours.]\n');

% % ========================================================
% % 2) Full cross-validation to find:
% %  - optimal state dimensionality for all methods
% %  - optimal smoothing kernel width for two-stage methods
% % ========================================================
% 
% % Select number of cross-validation folds
% numFolds = 10;
% 
% % Perform cross-validation for different state dimensionalities.
% % Results are saved in mat_results/runXXX/, where XXX is runIdx.
% for xDim = [2 5 8]
%   neuralTraj(runIdx, dat, 'method',  'pca', 'xDim', xDim, 'numFolds', numFolds);
%   neuralTraj(runIdx, dat, 'method', 'ppca', 'xDim', xDim, 'numFolds', numFolds);
%   neuralTraj(runIdx, dat, 'method',   'fa', 'xDim', xDim, 'numFolds', numFolds);
%   neuralTraj(runIdx, dat, 'method', 'gpfa', 'xDim', xDim, 'numFolds', numFolds);
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
% plotPredErrorVsDim(runIdx, kernSD);
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
% plotPredErrorVsKernSD(runIdx, xDim);
% % NOTES:
% % - This figure is used to find the optimal smoothing kernel for the
% %   two-stage methods.  The same smoothing kernel is used for all units.
% % - In this particular example, the optimal standard deviation of a
% %   Gaussian smoothing kernel with FA is 30 ms.
% % - Analogous to Figures 5B and 5C in Yu et al., J Neurophysiol, 2009.