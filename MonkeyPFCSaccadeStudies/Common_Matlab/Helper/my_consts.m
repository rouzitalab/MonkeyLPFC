function my_consts
%my_consts
%A script to load constants.

% global TOC              % structure: trialOutcome codes. Come from ptbmat
% global anaKeyword       % cellStructure: valid Keywords, number of Parameters
% global trialEvent       % structure: every possible Event has a special number and a textfield that can be evaluated

%% Natural constants
global CONST
CONST.dist1cmPerDeg = 57;

%% DEFine experimental/data parameter and event codes.
global DEF
% 
DEF.samplesPerSec         = 1000; % i.e. 
DEF.spikeTimeResolution   = DEF.samplesPerSec ; % Hz
DEF.spikeTimeFormat       = 'int32';
DEF.spikeSDFFormat        = 'single';
DEF.spikePatternFormat    = 'logical'; % logicals cannot be defined'in-line': use 'int8' for initialization and cast to 'logical'. That keeps the memory constant 
DEF.trialSpikeTimeFormat  = 'int32';  % change lim.maxTrialDuration accordingly!
DEF.contNeuroDataFormat   = 'single';
DEF.nsmatDataFormat       = 'int16';
DEF.eyeSampleRate         = 1000; %Hz never use this value directly, but transfer it to the global variable to allow an easy change if we ever have files with a different sampling rate
DEF.eyePosTimeFormat      = DEF.spikeTimeFormat;
DEF.trialEyePosTimeFormat = DEF.trialSpikeTimeFormat;
DEF.eyePosGazeFormat      = 'single';
DEF.eyePosPupilFormat     = 'int16';
DEF.eyePosDegreeFormat    = DEF.eyePosGazeFormat;
DEF.eyePosVelocityFormat  = DEF.eyePosGazeFormat;
DEF.eyePosResolutionFormat= DEF.eyePosGazeFormat;
DEF.eyePosTrialCountFormat= 'uint16';
DEF.eventTrialTimeFormat  = DEF.trialSpikeTimeFormat;
DEF.errorcodes            = 'int8';
% EYESTATEs
DEF.EYE_State_blink     = 0;
DEF.EYE_State_invalid   = 1;
DEF.EYE_State_offscreen = 2; 
DEF.EYE_State_undefined = 3;   
DEF.EYE_State_fixating  = 4; 
DEF.EYE_State_pusuit    = 5;
DEF.EYE_State_saccade   = 6;
% Invalidate Units
DEF.unitState.default         = -1;
DEF.unitState.valid           = 0; 
DEF.unitState.highBaseline    = 1;
DEF.unitState.manualExclusion = 2;
% anaKeywords 
DEF.anaKeyword(1).type = 'anaPeriod';               DEF.anaKeyword(1).numPars = 10;
DEF.anaKeyword(2).type = 'filterContNeuroData';     DEF.anaKeyword(2).numPars = 5;
DEF.anaKeyword(3).type = 'spikeTriggeredLFP';       DEF.anaKeyword(3).numPars = 8;
% trialEvents
DEF.trialEvent(1).type = 'trialStart';              DEF.trialEvent(1).eval = 'cDat.trial(t).startTime';
DEF.trialEvent(2).type = 'leverDown';               DEF.trialEvent(2).eval = 'cDat.trial(t).leverDown';
DEF.trialEvent(3).type = 'flipScreen';              DEF.trialEvent(3).eval = 'cDat.trial(t).flipScreen';
DEF.trialEvent(4).type = 'leverRelease';            DEF.trialEvent(4).eval = 'cDat.trial(t).leverRelease';
DEF.trialEvent(5).type = 'targetSaccade';           DEF.trialEvent(5).eval = 'eData.trial(t).targetSaccadeTime';  % based on Florian's valid saccade crit
DEF.trialEvent(6).type = 'gazeEnteredTargetTime';   DEF.trialEvent(6).eval = 'eData.trial(t).gazeEnteredTargetTime';
DEF.trialEvent(7).type = 'stopTime';                DEF.trialEvent(7).eval = 'eData.trial(t).stopTime';
DEF.trialEvent(8).type = 'saccadeTime';             DEF.trialEvent(8).eval = 'eData.trial(t).saccadeTime'; 

% Map flipScreens to events.
DEF.flipEvent = {'fixationOnset' 'targetOnset' 'cueOnset' 'cueOffset' 'fixationOffset' 'blank'};

DEF.parallelCompConfig = '';
if ~verLessThan('matlab', '7.7.0.471') %% 2008b
    toolboxVersions = ver;
    for i = 1:length(toolboxVersions)
       if strcmp(toolboxVersions(i).Name,'Parallel Computing Toolbox')
           if ~verLessThan('distcomp', '4.0')
               DEF.parallelCompConfig = 'local'; % use this configuration, may be different for dfferent machines
           end
       end
    end
end

DEF.fig_res_dpi = 600;
DEF.fig_width_cm = [8.9 18];% JNS: [8.5 11.6 17.6]; JNP: [8.9 18];

%% TrialFlag Definitions: second order qualification code for each trial,
% representing some stage of e.g. learning, average performance, block of
% datacollection, etc
global TFLAG            % structure: trialFlags
TFLAG.default   = -1;
TFLAG.undecided = 0;
TFLAG.learning  = 1;
TFLAG.learned   = 2;
TFLAG.learningAndLearned  = 3;

%% Limits on data, beyond which data are suspect.
global lim              % structure: limiting variables
lim.maxNumUnitsPerTrode = 5; % unsorted + 2 sorted units
lim.neuralChannelOffset = 1; % first 'hardware' channel to take into account 
lim.maxNumNeuroChannels = 32; % as we currently only have 32 channels
lim.maxNumTrodes        = lim.maxNumNeuroChannels; % if for any reason one wants to analyze tetrodes or multitrode
%Analog Channels
lim.maxNumUnitsPerAnalogChannel = 2; %thresholded + 1unit
lim.analogChannelOffset         = 129; % first analog channel
lim.maxNumAnalogChannels        = 16; % as we currently only have 32 neuro channels
%Continuous Channels
lim.contChannelOffset      = 147; % first possible continuous channel (ch 145 and 146 contain 'digital inputs' and 'serial inputs')
lim.maxNumContChannels     = lim.maxNumNeuroChannels+lim.maxNumAnalogChannels;
lim.contChannelSampleRate  = [1000]; % Hz
lim.contChannelReadBuffer  = 1000000;% samples ( == 8MB@ double read)
lim.contChannelFilterTypes = {'butter'};% butterworth
lim.contChannelFilterPoleFactor = 2; % polenumber must be a multiple of 2
lim.maxContChannelBands    = 5; % maximum number of Frequency bands that can recently be handeled 
lim.maxETAwaveforms        = 5000; % limits the number of events in an event triggered average to avoid unnecessary memory consumption 
% Trial
lim.maxNumTrials         = 3500; % trial
lim.maxTrialDuration     = min(100000,intmax('int32')); % ms : must not be larger than DEF.trialSpikeTimeFormat
lim.maxAnalysisPeriodDuration = min(10000,intmax('int32'));  % added AJS Aug 2011  

% %%
% global DIO              % structure: words that define specific events, i.e. leverDown etc.
% % Port 'A' - bits
% DIO.rewardBit       = 1; 
% DIO.unusedBit2      = 2;
% DIO.unusedBit3      = 3;
% DIO.unusedBit4      = 4;
% DIO.unusedBit5      = 5;
% DIO.unusedBit6      = 6;
% DIO.unusedBit7      = 7;
% DIO.strobeBit       = 8;
% % Port 'B' - Words
% DIO.startTrialWord      = 1;
% DIO.stopTrialWord       = 2;
% DIO.leverDown           = 3;
% DIO.leverRelease        = 4;
% DIO.flipScreen          = 5;
% DIO.startMyStim         = 6;
% DIO.fixating            = 7;
% DIO.breakFixation       = 8;
% 
% DIO.fixReqON            = 248;
% DIO.fixReqOFF           = 249;
% DIO.startSentenceWord   = 250;
% DIO.stopSentenceWord    = 251;
% DIO.startRecordingWord  = 252;
% DIO.stopRecordingWord   = 253;
% DIO.pauseRecordingWord  = 254;
% DIO.resumeRecordingWord = 255;
% 
% %% list of sentenceID's and their respective numbers
% global sentenceID       % structure: self explaining
% sentenceID.trialID = 1;
% sentenceID.fileID  = 2;

%% TrialOutcome Definitions should be delivered by .ptbmat-file !!!
% TOC.undef                = -1;
% TOC.hit                  = 0;
% TOC.noFixation           = 1;
% TOC.earlyFixationError   = 2;
% TOC.lateFixationError    = 3;
% TOC.earlyLeverRelease    = 4;
% TOC.anticipatedResponse  = 5;
% TOC.noResponse           = 6;
% TOC.notAttempted         = 7;
% TOC.notAttemptedLeverDown= 8;