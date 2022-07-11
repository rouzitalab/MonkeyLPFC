%Script to load, convert, consolidate, preprocess, then save.

%% Paths and Constants
addpath(genpath(fullfile(pwd, '..', '..', 'Common')));
my_consts;
my_paths;
sessions = my_sessions('TrickySessions');
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

%%
% JerryLee

% fmap sessions
%   ptb.eData.trial(ix)
%       .expType: 'M'
%       .attempt: some 0's, mostly 1's
%       .outcomeCode: 0:7
%       .targetChoice: 1:8; U,UR,R,DR,D,DL,L,UL
%                       9-16 are the same as 1-8 but greater distance
%       .trialoutcome: 0, 1, or empty. 1 if trial succeeded.
%       .flipScreen length 1:7
%       (newer) .cueVal: all 13's
%       (newer) .targetVal: all 21's
% Movement to targetChoice is made between flip 2 and flip 3.

% sra1 sessions
%   ptb.eData.trial(ix)
%       .expType: 'A'
%       .attempt: Almost all 1's, a few 0's
%       .outcomeCode: 0:12, mostly 0, 9, 12
%       .flipScreen len 1 to len 5
%       .targetChoice: 0's, 1's, even split. 0 is left, 1 is right.
%       .SRConfiguration : ones
%       (newer) .class : almost all 1,2, very few 13,14, 1=targ0, 2=targ1
%       Saccades typically between flip 3 (cueTime) and flip 4

% fmap sessions - mixed
%   ptb.eData.trial(ix)
%       .expType: 'A' or 'M'
%       .attempt: 0s and 1's
%       .outcomeCode: M: _0_,1,_2_,4,5,7; A: _0_,1,2,3,7,_9_,10,_12_
%       .flipScreen: M: 1,2,3,4,7; A: 1,2,3,4,5
%       .targetChoice: M: 0's, 1:8 even.; A: 0's 1's even.
%       .cueVal: M: all 13's; A: nothing
%       .targetVal: M: all 21's; A: nothing
%       .SRConfigures: M: nothing; A: ones
%       .class: M: nothing; A: 1,2 (even split)

% sra3_1 sessions
%   ptb.eData.trial(ix)
%       .expType: 'C'
%       .block: 1:etc.
%       .attempt: all 0's
%       .outcomeCode: 0 2 3 9 10 11 12; mostly 0, 3
%       .flipScreen: len 2 to len 6
%       .targetChoice: 1:8
%       .SR3InitialtrainingLevel: all 0's, all 3's.
%       .class: 41,45 12,16 26,30, 19,23 1,5 etc.

% sra3_2 sessions
%   ptb.eData.trial(ix)
%       .expType: 'D'
%       .attempt: all 0's
%       .outcomeCode: 0 2 3 9 10 11 12; mostly 0, some 3 and some 9
%       .flipScreen: len 2 to len 6
%       .targetChoice: all 2's in first file; 2 4 6 8 in second file. 1:8
%       in later files; empty in some files.
%       .class: 2 or 10, about even in first; 2, 6, 36, 40 (2/6 block 1, 36/40 block 2) in second file.
%               Lots of different classes with much higher numbers in later
%               files.
%       .class:.targetChoice ... 2:2; 6:6; 36:4; 40:8
%       (newer) .block: 1,2,etc.

% Marty
%
% sra3_1 sessions
%   ptb.eData.trial(ix)
%       .expType: 'C'
%       .block: 1:etc.
%       .attempt: all 0's
%       .outcomeCode: 0 2 3 9 10 11 12; mostly 0, 3
%       .flipScreen: len 2 to len 6
%       .targetChoice: 1:8
%       .SR3InitialtrainingLevel: all 0's, all 3's.
%       .class: 41,45 12,16 26,30, 19,23 1,5 etc. 1:46

% sra3_2 sessions
%   ptb.eData.trial(ix)
%       .expType: 'D'
%       .attempt: all 0's
%       .outcomeCode: 0 2 3 9 10 11 12; mostly 0, some 3 and some 9
%       .flipScreen: len 2 to len 6
%       .targetChoice: all 2's in first file; 2 4 6 8 in second file. 1:8
%                   : empty in late Marty files.
%       in later files; empty in some files.
%       .class: 2 or 10, about even in first; 2, 6, 36, 40 (2/6 block 1, 36/40 block 2) in second file.
%               Lots of different classes with much higher numbers in later
%               files.
%       .class:.targetChoice ... 2:2; 6:6; 36:4; 40:8
%       (newer) .block: 1,2,etc.