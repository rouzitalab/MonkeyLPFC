function [DIO, sentenceID] = get_default_DIO_sentenceID
DIO.startTrialWord      = 1;
DIO.stopTrialWord       = 2;
DIO.leverDown           = 3;
DIO.leverRelease        = 4;
DIO.flipScreen          = 5;
DIO.startMyStim         = 6;
DIO.fixating            = 7;
DIO.breakFixation       = 8;

DIO.fixReqON            = 248;
DIO.fixReqOFF           = 249;
DIO.startSentenceWord   = 250;
DIO.stopSentenceWord    = 251;
DIO.startRecordingWord  = 252;
DIO.stopRecordingWord   = 253;
DIO.pauseRecordingWord  = 254;
DIO.resumeRecordingWord = 255;
% list of sentenceID's and their respective numbers
sentenceID.trialID = 1;
sentenceID.fileID  = 2;