function output = procBehavBlocks(ptbTrials, useWindowTail, trialWindow, pResponded)
%output = procBehavBlocks(ptbTrials, useWindowTail, trialWindow, pResponded)
%Processes behaviour to return a table with information about behavioural
%performance. Each row in the table is a trial. Table columns are:
%TrialIndex, BlockIndex, Class, Target, Distractor, NewOCC, OCC, CueColor,
%Context, TruePositive, FalsePositive, FalseNegative, TrueNegative, IsGood
%Where TP, FP, FN, TN, IsGood are calculated using a trialWindow sliding
%window on the trials if:
%-all trials in the window belong to the same block (color-context pairing,
%should match [trial.block])
%-at least pResponded proportion of trials in block were completed
%-bias < maxBias (to make sure monkeys were not giving same response each
%trial).
if ~exist('useWindowTail','var') || isempty(useWindowTail)
    useWindowTail = true;
end
if ~exist('trialWindow', 'var') || isempty(trialWindow)
    trialWindow = 30;
end
if ~exist('pResponded', 'var') || isempty(pResponded)
    pResponded = 0.6;
end

%% Constants
minCorrPerStep = 2;
maxBias = 1.0;
pmarg = 0.01; %If probability approaches 0 or 1, fix by this amount.
colHeaders = {'TrialIndex', 'BlockIndex',...
    'Class', 'Target', 'Distractor',...
    'NewOCC', 'OCC', 'CueColor', 'Context',...
    'TrueA', 'FalseA', 'TrueB', 'FalseB', 'DPrime', 'BiasC', 'Acc', 'IsGood'};

%% Initialize table. TODO: Use Matlab table
nTrials = length(ptbTrials);
outputCell = cell(nTrials, length(colHeaders));

%% Shorthand some useful vectors.
tr_ix = [ptbTrials.expTrial]';
tr_blix = [ptbTrials.newBlock]';
tr_hiat = [ptbTrials.hiatus]';
tr_class = [ptbTrials.class]';
tr_targ = [ptbTrials.targClass]';
tr_dist = [ptbTrials.distClass]';
tr_newocc = [ptbTrials.newOutcomeCode]';
tr_occ = [ptbTrials.outcomeCode]';
tr_colours = {ptbTrials.cueColour}';
tr_context = {ptbTrials.targRule}';
tr_sac = [ptbTrials.sacClass]';

tr_ta = nan(nTrials, 1); tr_fa = tr_ta; tr_tb = tr_ta; tr_fb = tr_ta;
tr_c = tr_ta; tr_d = tr_ta; tr_acc = tr_ta; tr_isgood = tr_ta;

%% Use sliding window to calculate performance
nSteps = nTrials - trialWindow + 1;
tr_completed = (tr_newocc==0 | tr_newocc==9) & ~tr_hiat;
for ww = 1:nSteps
    winBool = false(nTrials,1); winBool(ww:ww+trialWindow-1)=true;
    isOneBlock = length(unique(tr_blix(winBool)))==1;
    isTrying = sum(tr_completed(winBool)) >= (pResponded*trialWindow);  % Monkey is trying to hit targets.
    if isOneBlock && isTrying
        bl = winBool & tr_completed;  % Completed trials in this block
        
        %Calculate ta,fa,tb,fb for this Window.
        tcs = unique(tr_context(bl));  % Should be exactly two
        
        % Can only calculate bias and d' if two options were available.
        if length(tcs) > 1
            targ_a = strcmpi(tr_context(bl), tcs{1});
            targ_b = strcmpi(tr_context(bl), tcs{2});  % Should be complementary to targ_a
            
            targ_corr = tr_newocc(bl) == 0;
            
            ta = sum(targ_a & targ_corr);
            fa = sum(targ_b & ~targ_corr);
            tb = sum(targ_b & targ_corr);
            fb = sum(targ_a & ~targ_corr);
            
            [d, c, acc] = performanceMeasures(ta, fa, fb, tb, pmarg);
            
            if useWindowTail
                w_ix = ww+trialWindow-1;  % = find(bl, 1, 'last')
            else
                w_ix = ww + round(trialWindow/2);
            end
            tr_ta(w_ix) = ta;
            tr_fa(w_ix) = fa;
            tr_tb(w_ix) = tb;
            tr_fb(w_ix) = fb;
            tr_d(w_ix) = d;
            tr_c(w_ix) = c;
            tr_acc(w_ix) = acc;
            % Mark trial as 'good' if ...
            % 1. abs(bias) is below maxBias AND...
            %    we have at least minCorrPerStep correct trials in each class
            % OR
            % 2. Only A trials, but none incorrect
            % 3. Only B trials, but none incorrect
            tr_isgood(w_ix) = ...
                (abs(c)<maxBias && ta>minCorrPerStep && tb>minCorrPerStep) ||...
                (fb == 0 && ta > 0) ||...
                (fa == 0 && tb > 0);
        end
        
    end
end

%% Bring vectors into output cell array
outputCell(:, 1:7) = mat2cell([tr_ix tr_blix tr_class tr_targ tr_dist tr_newocc tr_occ], ones(nTrials,1), ones(7,1));
outputCell(:, strcmpi(colHeaders, 'CueColor')) = tr_colours;
outputCell(:, strcmpi(colHeaders, 'Context')) = tr_context;
outputCell(:, 10:end) = mat2cell([tr_ta tr_fa tr_tb tr_fb tr_d tr_c tr_acc tr_isgood], ones(nTrials,1), ones(8,1));

output = cell2table(outputCell, 'VariableNames', colHeaders);