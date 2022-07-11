function [sacDir, sacAmp] = getSaccadeDirection(eData)
%sacDir = getSaccadeDirection(eData)
nTrials = length(eData.trial);
sacDir = nan(nTrials, 1);
sacAmp = sacDir;
for tt = 1:nTrials
    trial = eData.trial(tt);
    sacIx = 1;
    if ~isempty(trial.sacIx)
        sacIx = trial.sacIx;
    end
    sac = trial.saccades(sacIx);
    sacDir(tt) = sac.vectDir;
    sacAmp(tt) = sac.vectAmp;
end