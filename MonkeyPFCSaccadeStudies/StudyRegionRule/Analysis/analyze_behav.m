%Script to load preprocessed behavioural data then analyze that data.

%% Paths and Constants
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
global paths
sessions = my_sessions('RegionRule');
analysisParams = my_anaparams('behaviour');

sess_ix = 1;
%%
for sess_ix = 1:length(sessions)
    %% Load behavioural data.
    this_sess = sessions(sess_ix);
    pfnFull = fullfile(paths.preprocessed, 'ptb', this_sess.ptbfname);
    ptb = load(pfnFull, '-mat');

%     %% Plot the saccades for each trial
%     for tix = 1:length(ptb.eData.trial)  % 174 is a good example.
%         previewGazePerTrial(ptb.eData.trial(tix))
%         w = waitforbuttonpress;
%     end

    %% Re-evaluate behaviour.

    % Get the stimulus information for each trial, including the target
    % location.
    [ptb.eData.trial, ptb.params] = getTrialStimInfo(ptb.eData.trial, ptb.params);

    % Examine eye tracker data.
    % Reclassify saccades depending only on eye tracker data into
    % trial.newOutcomeCode: 0 = to target; 9 = to distractor; -1 for else.
    % Also updates each trial's saccade start time (sec) in
    % trial.sac_start_time.
    ptb.eData.trial = getTrialBehavResult(ptb.eData.trial, ptb.params, analysisParams);

    % Re-classify trials
    ptb.eData.trial = getNewClass(ptb.eData.trial, analysisParams);

    % Remove trials that are not of interest.
    % Note that, when we have cDat (merged with neural dat), we would use
    % triageTrials.m instead of the following line.
    ptb.eData.trial = ptb.eData.trial(ismember([ptb.eData.trial.newOutcomeCode], analysisParams.outcomeCodes));


    %% Plot saccade map.
    figure('Name', [this_sess.subject ' ' this_sess.ptbfname]);
    %plotSaccadeMap(ptb, analysisParams);
    
    behavOutput = procBehavBlocks(ptb.eData);
    plotBehav(behavOutput)

    %%
end