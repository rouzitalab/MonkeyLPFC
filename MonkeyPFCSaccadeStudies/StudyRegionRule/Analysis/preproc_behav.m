%Script to load, convert, consolidate, preprocess, then save.

%% Paths and Constants
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
sessions = my_sessions('RegionRule');
subjects = my_subs();
global paths

%% Step through each subject*session
for sess_ix = 1:length(sessions)
    this_sess = sessions(sess_ix);
    this_sub = subjects(strcmpi({subjects.name}, this_sess.subject));
    %% Load Pyschophysics Data
    ptb = getPTB(fullfile(paths.ptbedfRoot, this_sub.sname, 'ptbmat', this_sess.ptbfname));
    %% Load Eyetracker Data
    gazeData = getEyeData(this_sess.edfname, fullfile(paths.ptbedfRoot, this_sub.sname, 'edf'));
    %edfdata = edfmex(fullfile(paths.ptbedfRoot, name_map.(sub_names{sub_ix}), 'edf', efn));
    %1st converts from EDF if necessary
    %gazeData is a struct with fileds sample, gaze, pupil
    %The data are not yet divided into trials.

    %% Pre-process eyetracker data
    if isfield(ptb.eData, 'eyePosCalibAdjust')
        adjustments = ptb.eData.eyePosCalibAdjust; %TODO: account for differences in sampling rate.
    else
        adjustments = [];
    end
    gazeData = preProcessEyeData(gazeData, 'adjustments', adjustments, 'ptbParams', ptb.params);
    %Corrects for adjustments, smooths with gazeFilter,
    %converts .gaze (pixels) to .degree,
    %gets .state (keep large blocks of blink-free and on-screen)
    %.velocity, .acceleration, .dir, .velNorm
    clear adjustments

    %% Consolidate behavioural data (ptb.eData) and eye tracker data (gazeData)
    ptb.eData = consolidateBehavAndGaze(ptb.eData, gazeData);
    %Replaces eData.trial(n).eyePosition with the better data from gazeData.
    %Scans the eye tracker data from the whole trial (+small buffer) for saccades.

    %% Save this preprocesed data.
    if ~isdir(fullfile(paths.preprocessed, 'ptb'))
        mkdir(fullfile(paths.preprocessed, 'ptb'));
    end
    saveFull = fullfile(paths.preprocessed, 'ptb', this_sess.ptbfname);
    save(saveFull, '-struct', 'ptb');

    % Delete the ASC
    ascPath = fullfile(paths.preprocessed, 'asc');
    ascFilename = [this_sess.edfname(1:end-3) 'asc'];
    delete(fullfile(ascPath, ascFilename));

end
clear sess_ix this_sess this_sub ptb gazeData adjustments saveFull ascPath ascFilename

