function [trialBool, sacBool, tleBool, occBool, countBool] = triageTrials(cDat, varargin)
%cDat = triageTrials(cDat, varargin)
%varargin is name-value pairs
%'timeLockEvent' 'outcomeCodes' 'minTrialsPerGroup'
%Marks cDat.trial.invalid

%Default params
params = varg2params(varargin,...
    struct( 'anaWins', struct('timeLockEvent', 'saccadeOnset'),...
            'outcomeCodes', [0 9],...
            'dprimeLims', [-Inf Inf],...
            'minTrialsPerGroup', 20,...
            'distractorsDesired', 'any',...
            'minSacLatency', 200),...  % After cueOffset.
    {'anaWins' 'outcomeCodes' 'dprimeLims' 'minTrialsPerGroup', 'distractorsDesired' 'minSacLatency'});

%% Only grab trials with saccades
%otherwise, monkey attention/intention/processing undetermined.
sacBool = ~isnan([cDat.trial.sacStartTime]);

%% Limit to trials with all desired timeLockEvents
tles = unique({params.anaWins.timeLockEvent});
if any(ismember(tles, cDat.flipNames))
    [~, flip_ix] = ismember(tles, cDat.flipNames);
    [flip_ix, ix] = max(flip_ix);
    
    %Could check to make sure the trial had data past flip to accommmodate
    %all anawins... but nah.
    anawins = params.anaWins(strcmpi({params.anaWins.timeLockEvent}, tles{ix}));
    past_flip = max([anawins.winEdges]);
    
    nFlips = zeros(size(cDat.trial));
    for tr_ix = 1:length(cDat.trial)
        nFlips(tr_ix) = length(cDat.trial(tr_ix).flipScreen);
    end
    tleBool = nFlips >= flip_ix;
else
    tleBool = true(size(cDat.trial));
end
clear nFlips tr_ix

%% Only trials with appropriate newOutcomeCodes
if ~isempty(params.outcomeCodes)
    occBool = ismember([cDat.trial.newOutcomeCode], params.outcomeCodes);
else
    occBool = true(size(cDat.trial));
end

%% Only trials that were not on hiatus
hiBool = ~[cDat.trial.hiatus];

%% Remove trials with saccades too soon after cue extinguish
% Saccades are only supposed to come after fixation extinguish, but often
% the monkeys anticipated that and went too early. If we eliminate all the
% trials where saccades were earlier than 50 msec after fixation extinguish
% then we eliminate a lot of trials. To keep our trial count high, we
% eliminate trials where saccades were within 200 msec after cue
% extinguish (< 1% of trials).
saclatBool = true(size(sacBool));
if isfield(cDat, 'flipNames')
    flip_ix = find(strcmpi(cDat.flipNames, 'cueOffset'));
    if ~isempty(flip_ix)
        sac_trials = cDat.trial(sacBool);
        temp_bool = true(size(sac_trials));
        nsac = length(sac_trials);
        for tr_ix = 1:nsac
            temp_bool(tr_ix) = (sac_trials(tr_ix).sacStartTime - sac_trials(tr_ix).flipScreen(flip_ix)) > params.minSacLatency;
        end
        saclatBool(sacBool) = temp_bool;
    end
    clear flip_ix sasc_trials temp_bool nsac tr_ix
end

%% Remove trials that didn't mean DPrime criteria.
if params.dprimeLims(1) ~= -Inf
    dprimeBool = [cDat.trial.DPrime] >= params.dprimeLims(1);
else
    dprimeBool = true(size(cDat.trial));
end
if params.dprimeLims(2) ~= Inf
    dprimeBool = dprimeBool & [cDat.trial.DPrime] <= params.dprimeLims(2);
end

%% Remove trials in classes that have less than params.minTrialsPerGroup
%after eliminating trials thus far
uq_types = unique({cDat.trial.newType});
uq_classes = unique([cDat.trial.newClass]); %For this experiment
countBool = true(size(cDat.trial));
for ty_ix = 1:length(uq_types)
    ty_bool = strcmpi({cDat.trial.newType}, uq_types{ty_ix});
    for cl_ix = 1:length(uq_classes)
        cl_bool = [cDat.trial.newClass] == uq_classes(cl_ix);
        if sum(sacBool & tleBool & occBool & hiBool & ty_bool & cl_bool & saclatBool & dprimeBool) < params.minTrialsPerGroup
            countBool(ty_bool & cl_bool) = false;
        end
    end
end

%% Remove trials from a block that was really short
% or a block that did (not) have distractors inappropriately
if isfield(cDat.trial, 'cueTargRuleGroup')
    [count, group_ix] = hist([cDat.trial.cueTargRuleGroup], 1:max([cDat.trial.cueTargRuleGroup]));
    bad_block = group_ix(count < params.minTrialsPerGroup);
    blockBool = ~ismember([cDat.trial.cueTargRuleGroup], bad_block);
    
    if ~strcmpi(params.distractorsDesired, 'any')
        blk_distr_bool = nan(size(cDat.trial));
        for bl_ix = 1:length(group_ix)
            temp_bool = [cDat.trial.cueTargRuleGroup] == group_ix(bl_ix);
            blk_distr_bool(temp_bool) = any([cDat.trial(temp_bool).hasDistr]);
        end
        if strcmpi(params.distractorsDesired, 'no')
            blk_distr_bool = ~blk_distr_bool;
        end
        blockBool = blockBool & blk_distr_bool;
    end
else
    blockBool = true(size(sacBool));
end



%%
boolTable = [sacBool; tleBool; occBool; hiBool; countBool; blockBool; saclatBool; dprimeBool];
trialBool = sacBool & tleBool & occBool & hiBool & countBool & blockBool & saclatBool & dprimeBool;
fprintf(strcat('Of %i trials, removed %i total trials due to\n',...
    '-no saccades: %i (%i unique),\n',...
    '-truncated trial: %i (%i),\n',...
    '-missed target: %i (%i),\n',...
    '-hiatus: %i (%i),\n',...
    '-group inadequate: %i (%i),\n',...
    '-block too short: %i (%i)\n',...
    '-fast saccades: %i (%i)\n',...
    '-failed dprimeLims: %i (%i)\n'),...
    length(cDat.trial), sum(~trialBool),...
    sum(~sacBool), sum(all(boolTable(2:end, :)) & ~sacBool),...
    sum(~tleBool), sum(all(boolTable([1 3:end], :)) & ~tleBool),...
    sum(~occBool), sum(all(boolTable([1:2 4:end], :)) & ~occBool),...
    sum(~hiBool), sum(all(boolTable([1:3 5:end], :)) & ~hiBool),...
    sum(~countBool), sum(all(boolTable([1:4 6:end], :)) & ~countBool),...
    sum(~blockBool), sum(all(boolTable([1:5 7:end], :)) & ~blockBool),...
    sum(~saclatBool), sum(all(boolTable([1:6 8:end], :)) & ~saclatBool),...
    sum(~dprimeBool), sum(all(boolTable([1:7 9:end], :)) & ~dprimeBool));