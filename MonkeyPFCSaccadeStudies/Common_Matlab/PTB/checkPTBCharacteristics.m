%%
addpath(genpath(pwd));
my_consts; %global CONST DEF TFLAG lim
my_paths; %global paths
global paths

%%
sessions = getSessions();
sessions = sessions([sessions.isgood]);
[tmp ind] = sort([sessions.date]);
temp_sessions = sessions(ind);
name_map.JerryLee = 'j'; %For the EDF and PTB files, a different name was used.
name_map.Marty = 'm';
fieldNames = {};
eFieldNames = {};
trialFieldNames = {};
paramFieldNames = {};

%Go through each file to find the unique field names and n_rows
fprintf('Finding all field names ');
for row_ix = 1:length(temp_sessions)
    ptb = load(fullfile(paths.ptbedfRoot, name_map.(temp_sessions(row_ix).subject), 'ptbmat', temp_sessions(row_ix).ptbfname), '-mat');
    n_rows = n_rows + 1;
    fieldNames = union(fieldNames, fieldnames(ptb));
    eFieldNames = union(eFieldNames, fieldnames(ptb.eData));
    trialFieldNames = union(trialFieldNames, fieldnames(ptb.eData.trial(1)));
    paramFieldNames = union(paramFieldNames, fieldnames(ptb.params));
    if mod(row_ix, 5)==0
        fprintf('.');
    end
end
fprintf(' Done.\n');

% Go through the files to identify which files have which fields
n_rows = length(temp_sessions);
expCodes = cell(n_rows, 1);
fieldPresent = false(n_rows, length(fieldNames));
eFieldPresent = false(n_rows, length(eFieldNames));
trialFieldPresent = false(n_rows, length(trialFieldNames));
paramFieldPresent = false(n_rows, length(paramFieldNames));
fprintf('Populating fieldPresent, eFieldPresent, trialFieldPresent, and paramFieldPresent. ');

%Hijacking this to get values of a couple fields
[temp_sessions.SR3errorStrategy] = deal(nan);
[temp_sessions.SR3trainingLevel] = deal(nan);
[temp_sessions.SR3InitialtrainingLevel] = deal(nan);

for row_ix = 1:length(temp_sessions)
    ptb = load(fullfile(paths.ptbedfRoot, name_map.(temp_sessions(row_ix).subject), 'ptbmat', temp_sessions(row_ix).ptbfname), '-mat');
    expCodes{row_ix} = ptb.params.expCode;
    fieldPresent(row_ix, :) = ismember(fieldNames, fieldnames(ptb));
    eFieldPresent(row_ix, :) = ismember(eFieldNames, fieldnames(ptb.eData));
    trialFieldPresent(row_ix, :) = ismember(trialFieldNames, fieldnames(ptb.eData.trial(1)));
    paramFieldPresent(row_ix, :) = ismember(paramFieldNames, fieldnames(ptb.params));
    if mod(row_ix, 5)==0
        fprintf('.');
    end
    
    %Hijack
    if isfield(ptb.params, 'SR3errorStrategy')
        temp_sessions(row_ix).SR3errorStrategy = ptb.params.SR3errorStrategy;
    end
    if isfield(ptb.params, 'SR3trainingLevel')
        temp_sessions(row_ix).SR3trainingLevel = ptb.params.SR3trainingLevel;
    end
    if isfield(ptb.eData.trial, 'SR3InitialtrainingLevel')
        temp_sessions(row_ix).SR3InitialtrainingLevel = unique([ptb.eData.trial.SR3InitialtrainingLevel]);
    end
end
fprintf(' Done.\n');

%% ptb
% First 6 files Have CALADJ, responseF, responseFM, responseR, responseS but not TOC
% The rest have TOC.
% All other fields (DIO, eData, params, sentenceID) are common.

%% ptb.eData
% startEyeRecordings -> firstEyeSample from file 4 on.
% -eyePosCalibAdjust is added 6:end
% -regularExpEnd and trialAlignentEvent are added 7:end
% -CLUT is added 8:end
% -newedfFile is added 59:end

%% ptb.eData.trial(1)
temp = [zeros(1, size(trialFieldPresent, 2)); diff(trialFieldPresent)];
for i = 1:size(temp,1)
    fprintf('\nRow %i:\n', i);
    col_bool = abs(temp(i,:))>0;
    if any(col_bool)
        fprintf('%s\n', datestr(temp_sessions(i).date));
        col_id = find(col_bool);
        for cc = 1:length(col_id)
            fprintf('%i %s\n', temp(i, col_id(cc)), trialFieldNames{col_id(cc)});
        end
    end
end

% -expTrial, expType, words are always present

% 19-Nov-2008:
% -response
% +TargetTime, cueTime, distractor2Time, eyeSyncStartTime, eyeSyncStopTime,
% eyeposition, flipScreen, outcomeCode, reason, targetChoice,
% timeleverdown, timeleverup, timeofrelease, trialattempt, trialoutcome

% 20-Nov-2008
% eyeSyncStartTime -> eyeSyncTime
% trialStartTime -> startTime
% trialStopTime -> stopTime
% trialID -> ID
% -eyeSyncStopTime

% 08-Dec-2008
% -TargetTime, timeleverdown, timeleverup
% +attempt, cueDelay, cueVal, distractorDelay, eyePosition, leverDownTime,
% leverReleaseTime, targetDelay, targetTime, targetVal

% 09-Dec-2008
% timeofrelease -> leverRelease
% trialattempt, trialoutcome

% 24-Dec-2008
% -cueDelay, cueVal, distractor2Time, distractorDelay, leverReleaseTime,
% reason, targetDelay, targetTime, targetVal
% +cuePresentTime, fix2cueDelay, gazeWithinTargetTime, rewardTime,
% time2HitTarget, time2StayInTarget

% 05-Jan-2009
% +SRConfiguration

% 13-Jan-2009
% cueTime -> cuePresentedTime
% +class, cueExtinguishedTime

% 29-Apr-2009
% +block

% 17-Jun-2009
% +SR3InitialtraininLevel

% From 08-May-2009, SR3errorStrategy is present in only SOME files.
% JL sra3_2 (task D) 2009-May, 08, 13, 20, 20, 22, 22-Jun
% M sra3_1 (task C) 16-Jul-2009 to 28-Sep-2009; sra3_2 (task D) 04-Oct-2009
% to 06-Nov-2009

% From 02-Jul-2009 SR3trainingLevel is present in only SOME files.

%% ptb.params
temp = [zeros(1, size(paramFieldPresent, 2)); diff(paramFieldPresent)];
for i = 1:size(temp,1)
    fprintf('\nRow %i:\n', i);
    col_bool = abs(temp(i,:))>0;
    if any(col_bool)
        fprintf('%s\n', datestr(temp_sessions(i).date));
        col_id = find(col_bool);
        for cc = 1:length(col_id)
            fprintf('%i %s\n', temp(i, col_id(cc)), paramFieldNames{col_id(cc)});
        end
    end
end
% Always present:
% FMExp, FeatureExp, ISI, Interleave, RadialExp, SpatialExp,
% analysedTrialTypes, anticipatedResponseDuration, backgroundIndexColor,
% basicDataSavePath, beep, button2cueTime, calibration, collectSpikes,
% cuePresentTime, currElectrodeConfig, currSession, delay2fixation,
% earlyResponsePermitted, expCode, expInfo, eyeAverage, filename,
% finalDataSavePath, fixChangeluminance, fixRad, giveFeedback, hitDelay,
% keyboardResponse, leverDownRewardDuration, leverReleaseDuration, lowbeep,
% missDelay, nCa, numberofTrials, preLeverBeepDelay,
% preLeverProgressDisplay, provideWords, refreshRate, responseTimeWindow,
% reward, rewardDuration, rewardbeep, sampleSource, saveData,
% spikeRecSystem, spikeSource, stimulator, subject, tempEDFName, time,
% usedFileNamesFile, waitForLever

% 19-Nov-2008
% +subjectScreenDistance, subjectScreenPixelPerDegree,
% subjectScreenResolution, subjectScreenSize, unitSubjectScreenDistance

% 08-Dec-2008
% +FMnumStimPerAnnulus

% 24-Dec-2008
% FM -> SR1
% -CueDelayDifficulty, DistractorDifficulty, FMCueDifficulty,
% FMTargetDifficulty, FMbutton2DistractorTime, FMcue2TargetTime,
% FMcuePresentTime, FMdistractor2CueTime, FMgridLength, FMnumAnnuli,
% FMnumStimPerAnnulus, FMptSize, TargetDelayDifficulty, Xtranslate,
% Ytranslate, fixationRequired
% +SR1Exp, SR1Xtranslate, SR1Ytranslate, SR1cueTime, SR1fix2cueDelay,
% SR1fixationRequired, SR1time2HitTarget, SR1time2StayInTarget, blue,
% cueHeight, extralowbeep, fixPointSize, fixPointfixRadius, green,
% penWidth, probUncertainTarget, red, targetDistance, targetHeight,
% targetPointfixRadius, yellow

% 05-Jan-2009
% +analysisBlockSize

% 22-Jan-2009
% +Add back some for FM
% 03-Feb-2009
% Remove some for FM

% 27-Apr-2009
% SR1 -> SR3
% +SR2Exp, SR3_1_Exp, SR3_2_Exp,
% black, white, cueRow, distractorColor, fixPointColor, initFixColor, targetColor