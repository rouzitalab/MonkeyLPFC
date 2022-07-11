function dat = as_SimSaccadeExp(nTrials,trialDur);
% dirVec = [1 2 3 4 5 6 7 8];

for iTrials=1:nTrials
    randDir = ceil(rand*8);
    data = as_simEnsembleSaccade(randDir,trialDur);
    spikeTrain = zeros(1,trialDur*1000);
    dat(iTrials).trialId = iTrials;
    spikes = zeros(32,trialDur*1000);
    for iCell = 1:32
        spikes ( iCell,data.spikeTimes{iCell}*1000 ) = 1;
    end
    dat(iTrials).spikes=spikes;
    dat(iTrials).dir = randDir;
    clear spikes
end