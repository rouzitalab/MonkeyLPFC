clear
clc
path = 'G:\Projects\gpfa_v0203\data\chadOutput\sra3_1_m_074_00+01_1_gpfa.mat';
data1 = load(path);

seqTrain = data1.trialTraj;
y = data1.trialCode;
col = [[1 0 0]; [0 1 0]; [0 0 1]; [0 1 1]; [1 1 0]; [1 0 1];[0.5 0.1 0.1];
    [0.1 0.5 0.1]; [0.5 0.5 0.5]; [0.5 0.5 1];[1 0.5 0.5]; [0.3 0.2 0.4]];
figure;
for n = 1:length(seqTrain)
    dat = seqTrain(n, :, [1, 2]);
    plot(dat(:,1), dat(:,2), 'o','color', col(y(n),:));
    hold on;
end
