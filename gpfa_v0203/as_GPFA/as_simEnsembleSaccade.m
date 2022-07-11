function data = as_simEnsembleSaccade(targetDir)

% currently assumes 1 second epoch with 1000 Hz
'test'
BWVar = 5000;
epochDur = 1000; % 1000 msec
sampleRate = 1000;
data.targetDir=targetDir;
if targetDir == 8
    targetDir == -2;
elseif targetDir == 7
    targetDir == -3;
elseif targetDir == 6;
    targetDir == -4;
end

ensembleSize=32;
ensembleTuningDir = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7 8 8 8 8];

if targetDir == 1
    distArray = [0 1 2 3 4 3 2 1];
elseif targetDir == 2
    distArray = [1 0 1 2 3 4 3 2];
elseif targetDir == 3
    distArray = [2 1 0 1 2 3 4 3];
elseif targetDir == 4
     distArray = [3 2 1 0 1 2 3 4];
elseif targetDir == 5
     distArray = [4 3 2 1 0 1 2 3];
elseif targetDir == 6
     distArray = [3 4 3 2 1 0 1 2];  
elseif targetDir == 7
     distArray = [2 3 4 3 2 1 0 1];
elseif targetDir == 8
     distArray = [1 2 3 4 3 2 1 0];
end
     

% elseif targetDir ==2;
%     ensembleTuningDir = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 -4 -4 -4 -4 -3 -3 -3 -3  -2 -2 -2 -2];
% elseif targetDir == 3;
% ensembleTuningDir = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 -4 -4 -4 -4 -3 -3 -3 -3  -2 -2 -2 -2];

ensembleTuningBandwidth = BWVar*[0.1 0.5 1 2 0.1 0.5 1 2 0.1 0.5 1 2 0.1 0.5 1 2 0.1 0.5 1 2 0.1 0.5 1 2 0.1 0.5 1 2 0.1 0.5 1 2 ];
ensembleLambda = 100*(ones(1,32));

for i=1:ensembleSize
    %     dist = abs(targetDir-abs(ensembleTuningDir(i))) ;
    %     if dist > 4
    %         dist = abs(4-abs(targetDir-abs(ensembleTuningDir(i))));
    %     end
    dist = distArray(ensembleTuningDir(i));
    k = normpdf(dist,0,ensembleTuningBandwidth(i))* 1/normpdf(0,0,ensembleTuningBandwidth(i));
    lambda = ensembleLambda(i)/sampleRate; % since 1000 msec in 1 sec
    spikeTrain = poissrnd(k*lambda,1,epochDur);
    spikeTimes = find(spikeTrain)/sampleRate;
    data.spikeTimes{i} = spikeTimes;
    data.numSpikes{i} = length(spikeTimes);
    data.distanceFromPreferrred{i} = dist;
end
    