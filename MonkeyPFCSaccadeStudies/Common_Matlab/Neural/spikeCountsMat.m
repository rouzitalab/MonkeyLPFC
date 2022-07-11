function [X, tVec] = spikeCountsMat(cDat, varargin)
%TODO: FINISH
%% Get time vector
% Currently does not allow trials that only participate in one window.

%Identify time limits for spike count bins
%starts = [cDat.trial.startTime] - [cDat.trial.alignTime] - [cDat.trial.eventTime];
%stops = [cDat.trial.stopTime] - [cDat.trial.alignTime] - [cDat.trial.eventTime];
starts = -1*timeFac.*[cDat.trial.eventTime];  % Each trial starts at t=0, or eventTime msec before the timeLockEvent.
stops = timeFac.*[cDat.trial.stopTime] + starts;

trialEdges = [max([params.winEdges(1) starts]) min([params.winEdges(end) stops])];

%Time vector with steps of binDuration
if trialEdges(1)<0
    t1 = 0:params.binWidth:abs(trialEdges(1));
    t1 = fliplr(-t1);
else
    t1 = trialEdges(1);
end
t2 = max(t1)+params.binWidth:params.binWidth:trialEdges(2);
edgeVec = [t1 t2];
tVec = edgeVec(2:end) - params.binWidth/2; %Bin Centers
nTimes = length(tVec);
clear trialEdges t1 t2