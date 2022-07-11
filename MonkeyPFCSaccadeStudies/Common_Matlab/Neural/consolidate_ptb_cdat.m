function cDat = consolidate_ptb_cdat(ptb, cDat)
%cDat = consolidate_eData_cDat(ptb, cDat)
jitterTol = 4;  % msec
slopeTol = 0.01;  %

%% Match trials
%Most likely is that there is one more eData.trial than cDat.trial
%This happens when the PTB starts a trial but the user interrupts it.
while length(ptb.eData.trial) > length(cDat.trial)
    ptb.eData.trial(end) = [];
end
if length(ptb.eData.trial) ~= length(cDat.trial)
    error('Cannot match trial lengths.');
end

%% Copy over information from eData to cDat
% eData.trial(t).flipScreen vs find(eData.trial(t).sentWords(1.:) == DIO.flipScreen)
% - temporal offset between events (should be ~ 4ms , check especially LEVER!!! DOWN and UP!!

X = num2cell([ptb.eData.trial.newClass]);
[cDat.trial.newClass] = X{:};

X = {ptb.eData.trial.newClassStr};
[cDat.trial.newClassStr] = X{:};

X = num2cell([ptb.eData.trial.newOutcomeCode]);
[cDat.trial.newOutcomeCode] = X{:};

X = num2cell(1000*[ptb.eData.trial.sacStartTime]);
[cDat.trial.sacStartTime] = X{:};  % Time in msec since start of trial.

X = num2cell(cat(1, ptb.eData.trial.sacPol), 2);
[cDat.trial.sacPol] = X{:};

X = num2cell(cat(1, ptb.eData.trial.targPol), 2);
[cDat.trial.targPol] = X{:};

X = num2cell(cat(1, ptb.eData.trial.sacEndXY), 2);
[cDat.trial.sacEndXY] = X{:};

X = {ptb.eData.trial.newType};
[cDat.trial.newType] = X{:};

X = num2cell([ptb.eData.trial.newBlock]);
[cDat.trial.cueTargRuleGroup] = X{:};

X = num2cell([ptb.eData.trial.hiatus]);
[cDat.trial.hiatus] = X{:};

X = num2cell(~cellfun(@isempty, {ptb.eData.trial.distStr}));
[cDat.trial.hasDistr] = X{:};

if isfield(ptb.eData.trial, 'DPrime')
    X = num2cell([ptb.eData.trial.DPrime]);
    [cDat.trial.DPrime] = X{:};
end

cDat.flipNames = ptb.params.flipNames;

%% Check to see if each trial's flipScreen's roughly match up.
[cDat.trial.flipBool] = deal(true);
for tt = 1:length(cDat.trial)
    cflip = cDat.trial(tt).flipScreen;
    eflip = 1000*ptb.eData.trial(tt).flipScreen;
    if length(cflip)==length(eflip)
        X = [ones(size(cflip)) eflip'];
        B = X\cflip;
        slopeOK = abs(B(2)-1) <= slopeTol;
        jitterOK = ~any( abs(cflip - X*B) > jitterTol);
        cDat.trial(tt).flipBool = slopeOK && jitterOK;
    else
        cDat.trial(tt).flipBool = false;
    end
end