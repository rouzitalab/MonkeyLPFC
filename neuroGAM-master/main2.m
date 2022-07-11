%% Loading the data
clear;
clc;
path = 'G:\Projects\neuroGAM-master\data';
data_list = dir(fullfile(path, '*.mat'));
for sess = 2:9
clearvars -except 'sess' 'data_list' 'path'
fprintf('Working on Sess %d\n', sess);
data_path = append(path, '\',data_list(sess).name);
data1 = load(data_path);
%% Limiting data to learned or not learned trials
%(finding the indices to apply later)
% [learned,not_learned] = monkey_performance(data1);
learned = data1.learned + 1;
unlearned = data1.unlearned + 1;
flag = learned;
% flag = unlearned;
% the expected best model
% best_model = 3;
best_model = 1;
saccade = data1.saccade;
target = data1.target;
spikes = data1.spikes;
times = data1.times;
% idx = find (times < 1.451);
% spikes = spikes(:,1:1850,:); % limit spikes to pre-saccade activities (J)
% spikes = spikes(:,1:132,:); % limit spikes to pre-saccade activities (M)
% spikes = spikes(:,idx,:);
cue = data1.color;
outcome = data1.outcome;
% flag = find(outcome~=-1);
% tmp = outcome(flag);
% flag2 = find(tmp==0);
spikes(isnan(spikes))=0; % zero-ing the NaN spikes
%% Extracting saccade, color, rule, and the spikes
color = double.empty();
color(length(cue)) = 0;
for i=1:length(cue) % assigning labels to the cue color : Red = 1, Green = 2, Blue = 3
    if isequal(cue{i},'r')
        color(i)=1;
    elseif isequal(cue{i},'g')
        color(i)=2;
    elseif isequal(cue{i},'b')
        color(i)=3;
    end
end


% applying learned or not learned indices
color = color(flag);
% rule = rule(flag);
% pair = pair(flag);
saccade = saccade(flag) + 1;
target = target(flag) + 1;
spikes = spikes(flag,:,:);
% color = color(flag2);
% saccade = saccade(flag2);
% spikes = spikes(flag2,:,:);
% tmp = permute(spikes,[2,1,3]);
% tmp = downsample(tmp,10);
% spikes = permute(tmp, [2,1,3]);
spikes = abs(spikes);
spikes = uint16(spikes);
%% getting variables ready for the neuroGAM method
trials = size(spikes,1); % number of time trials
samples = size(spikes,2); % number of time samples
t_trg = find(times==0);
t_cue = find(times==0.25);
t_mov = find(times==1.25);
% channels = size(spikes,3); % number of time channels
sac = zeros(1, trials * samples); % variables are in 1D format
trg = zeros(1, trials * samples);
col = zeros(1, trials * samples);
rul = zeros(1, trials * samples);
for i = 1:trials
    sac((i-1)*samples+1:i*samples) = saccade(1,i);
    trg((i-1)*samples+1:i*samples) = target(1,i);
    col((i-1)*samples+1:i*samples) = color(1,i);
end
% combined rule labels using saccade (or target) and the color cue
% direction = sac;
direction = trg;
for i = 1:trials*samples
    if direction(i)==1 && col(i)==1, rul(i)=1;
    elseif direction(i)==1 && col(i)==2, rul(i)=2;
    elseif direction(i)==1 && col(i)==3, rul(i)=3;
    elseif direction(i)==2 && col(i)==1, rul(i)=4;
    elseif direction(i)==2 && col(i)==2, rul(i)=5;
    elseif direction(i)==2 && col(i)==3, rul(i)=6;
    elseif direction(i)==3 && col(i)==1, rul(i)=7;
    elseif direction(i)==3 && col(i)==2, rul(i)=8;
    elseif direction(i)==3 && col(i)==3, rul(i)=9;
    elseif direction(i)==4 && col(i)==1, rul(i)=10;
    elseif direction(i)==4 && col(i)==2, rul(i)=11;
    elseif direction(i)==4 && col(i)==3, rul(i)=12;
    elseif direction(i)==5 && col(i)==1, rul(i)=13;
    elseif direction(i)==5 && col(i)==2, rul(i)=14;
    elseif direction(i)==5 && col(i)==3, rul(i)=15;
    elseif direction(i)==6 && col(i)==1, rul(i)=16;
    elseif direction(i)==6 && col(i)==2, rul(i)=17;
    elseif direction(i)==6 && col(i)==3, rul(i)=18;
    elseif direction(i)==7 && col(i)==1, rul(i)=19;
    elseif direction(i)==7 && col(i)==2, rul(i)=20;
    elseif direction(i)==7 && col(i)==3, rul(i)=21;
    elseif direction(i)==8 && col(i)==1, rul(i)=22;
    elseif direction(i)==8 && col(i)==2, rul(i)=23;
    else, rul(i)=24;
    end
end
tmp = unique(rul);
tmp2 = 1:length(tmp);
tmp3 = zeros(length(rul),1);
for i = 1 : length(rul)
    for j = 1 : length(tmp)
        if rul(i)==tmp(j), tmp3(i)=tmp2(j); end
    end
end
    
sac = sac';
trg = trg';
rul = tmp3;
col = col';

% the xt is the cell combining 1D labels
% SINCE IT IS NECESSARY TO HAVE THE SAME 1D DIMENSIONS FOR SPIKES AND
% VARIABLES: the labels are repeated #(time samples) * #(channels) times
% for each trial ~ 10,000 times for each trial.
xt = {sac; col; trg; rul};
% the prs structure for the method
prs.varname = {'Saccade','Color','Target','Rule'};
prs.vartype = {'1D','1D','1D', '1D'};
prs.basistype = {'raisedcosine','raisedcosine','raisedcosine','raisedcosine'};
prs.nbins = {8,3,8,length(unique(rul))}; % number of different labels in our case
prs.binrange = [];
prs.nfolds = 10;
prs.dt = 0.001; % 1/sampling_rate
prs.filtwidth = 10;
prs.linkfunc = 'log';
prs.lambda = {10, 10, 10, 10};
prs.alpha = 0.05;
prs.varchoose = [0 0 0 0]; % auto choice which variable is encoded by yt
prs.method = 'fastforward';
%% Calling the method on single channels
channels = size(spikes,3);
models = cell(channels,1);
best_channels = [];
for ch=1:channels
    fprintf(['Processing Channel ' int2str(ch) ' of ' int2str(channels) '\n']);
    spk = zeros(1, trials * samples);
    tmp = spikes(:,:,ch);
    for i=1:trials
        spk((i-1)*samples+1:i*samples) = tmp(i,:);
    end
    % the yt is the 1D spikes
    yt = spk';
    [models{ch}, ~, ~] = BuildGAM(xt,yt,prs);
    if models{ch}.bestmodel==best_model, best_channels(end+1)=ch; j = j+1; end
end
% %% Plotting all single channels
% for ch=1:channels
%     PlotGAM(models{ch}, prs);
% end
%% Plotting the best channels
% for i=1:length(best_channels)
%     j = best_channels(i);
%     PlotGAM(models{j}, prs);
% end
%% Free-up memory used by temporary variables
clear i k tmp tmp2 tmp3;
%% Write to excel file
printLog();
end
beep