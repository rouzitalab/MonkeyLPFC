function dat10 = as_shuffleTrajectories_dat1(dat)

nTrials= size(dat,2);
shuffleIdx = shuffle(1:nTrials);
for i=1:nTrials
    dat10(i).dir = dat(i).dir;
    dat10(i).trialId = dat(i).trialId;
    dat10(i).spikes = dat(shuffleIdx(i)).spikes;
end

