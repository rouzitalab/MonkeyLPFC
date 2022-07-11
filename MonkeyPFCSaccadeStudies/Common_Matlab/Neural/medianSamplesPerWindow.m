function [win_samples, win_length] = medianSamplesPerWindow(sessions, analysisParams, dhParams)
% Count how many samples per window per trial across sessions.
% This is needed so we know how much to timewarp each trial x window
global paths
win_length = cell(1, length(sessions));
for sess_ix = 1:length(sessions)
    %%
    this_sess = sessions(sess_ix);
    
    %% Load and process ptbmat
    ptb = load(fullfile(paths.preprocessed, 'ptb', this_sess.ptbfname), '-mat');
    [ptb.eData.trial, ptb.params] = getTrialStimInfo(ptb.eData.trial, ptb.params);
    ptb.eData.trial = getTrialBehavResult(ptb.eData.trial, ptb.params, analysisParams);
    ptb.eData.trial = getNewClass(ptb.eData.trial, analysisParams);
    %Add behavOutput.DPrime to ptb.eData.trial
    behavOutput = procBehavBlocks(ptb.eData.trial);
    DPrime = num2cell(behavOutput.DPrime);
    [ptb.eData.trial.DPrime] = DPrime{:};
    
    %% Load and process cerbDat
    cDat = load(fullfile(paths.preprocessed, 'cDat', [this_sess.edfname(1:end-3) 'mat']));
    cDat = consolidate_ptb_cdat(ptb, cDat);
    trialBool = triageTrials(cDat, analysisParams);
    cDat.trial = cDat.trial(trialBool);
    nTrials = length(cDat.trial);
    cDat = triageUnits(cDat, analysisParams);
    
    %% Count samples per trial
    win_length{sess_ix} = int16(zeros(length(analysisParams.anaWins), nTrials));
    for win_ix = 1:length(analysisParams.anaWins)
        ana_win = analysisParams.anaWins(win_ix);
        win_cDat = trimTrials(cDat, ana_win);
        win_D = cDat2DataHigh(win_cDat.trial);
        for tr_ix = 1:nTrials
            win_length{sess_ix}(win_ix, tr_ix) = ...
                floor(size(win_D(tr_ix).data, 2) / dhParams.binWidth);
        end
    end
end
win_samples = median(cat(2, win_length{:}), 2);