function [trials, ptbParams] = getTrialStimInfo(trials, ptbParams)
%function trials = getTrialStimInfo(trial, ptbParams)
%To each trial, adds fields
%.newType in {'M', 'A', 'SR3', 'SR4'}
%.cueColour in {'r', 'g', 'b'}
%(the following fields may be NaN if trial type does not have targ/dist)
%.targClass in 1:8 (or 1:16 for some 'M')
%.targXY = [xpos ypos] in pixels. Already includes centrePt offset.
%.targPol = [theta radius], from 0. Need to add centrePt offset after.
%.distClass
%.distXY
%.distPol

%TODO: For some reason, SR4 sometimes has SR3-like stimuli.
%TODO: Use data in flipScreen to create trial fields with common names.
%TODO: Handle cases when errorStrategy changes within a session.
%TODO: ExpType 'A'

%% Constants

%Cue Info. The PTB code uses trial classes for which the value indexes into
%a predefined colour-map, below.
cueColourMap = {'r', 'g', 'b'};
cueMatrix = [1 2; 1 3; 2 1; 2 3; 3 1; 3 2];  % 1=red; 2=green; 3=blue

%Target Info for SR3
%In Adam's code, the 1:8 part of the class is mapped to theta with
%targetTheta = 45*((1:8)-8-1); then the targetPt is calculated as
%target_yx = centre_yx + radius*[-cosd(targetTheta) sind(targetTheta)];
% Notice yx, not xy. Then he fliplr's the result.
%This can be more simply represented the following way.
targetTheta = deg2rad([270 315 0 45 90 135 180 225]);
targetStr = {'UU' 'UR' 'RR' 'DR' 'DD' 'DL' 'LL' 'UL'};  % The y-zero was up, so flip u/d.
centrePt = ptbParams.subjectScreenResolution / 2;

% Only needed for SR4: groups for general rules.
context_options = {'AnyUp' 'AnyDown' 'AnyLeft' 'AnyRight'};
% [d;d;d],[u;u;u]
% [l;l;l];[r;r;r]
targGroups = cat(3, [[8;1;2] [6;5;4]], [[8;7;6] [2;3;4]]);
targDistGroups = cat(4, targGroups, flip(targGroups, 2)); clear targs %

%%
for expType = unique({trials.expType})
    
    %% Prepare
    trialBool = strcmpi({trials.expType}, expType)';
    classes = [trials(trialBool).class]';
    nTrials = sum(trialBool);
    
    %% Get the trial type for all trials.
    trialTypes = repmat({'unknown'}, 1, nTrials);
    
    % Trial type can be determined based on the experiment type for expTypes M and A.
    % For expType C or D, we need more information, and that can change on a per-trial basis.
    
    if strcmpi(expType, 'M')
        radius = ptbParams.FMgridLength./(2*(ptbParams.FMnumAnnuli:-1:1));
        fixJitterXY = [0 0]; %Add to fix/targ/distr XY
        trialTypes = repmat({'M'}, 1, nTrials);
    elseif strcmpi(expType, 'A')
        radius = ptbParams.targetDistance;  % ??
        trialTypes = repmat({'A'}, 1, nTrials);
    elseif strcmpi(expType, 'C') || strcmpi(expType, 'D')
        % Some variables contain two indices, one for expType C and one for
        % expType D.
        tl_colix = 1;
        if strcmpi(expType, 'D')
            tl_colix = 2;
        end
        radius = ptbParams.SR3radius(tl_colix);
        fixJitterXY = [ptbParams.SR3Xtranslate(tl_colix) ptbParams.SR3Ytranslate(tl_colix)];
        
        % Trial type for expType C and D is more complicated.
        % Ultimately it comes down to whether the PTB code calls
        % dSR3T2c_setscreens2 or dSR4T2c_setscreens2. The code path is as
        % follows:
        % (start in saccadeGUImain)
        % if on probation -> doSR3Trial2e
        % if not on probation -> doSR3Trial2d
        %
        % Whether or not we are on probation can be changed throughout the
        % experiment with a key press -> saved to trials.SR3errorStrategy,
        % but only on the trial for which the key was pressed (how can we
        % know what the error strategy was at the beginning of the
        % experiment?!)
        %
        % doSR3Trial2e runs the same trial over and over until 2/3 are
        % successful, then it calls dSR3T2c_setscreens2
        % doSR3Trial2d runs dSR3T2c_setscreens2 if expType C, or
        % dSR4T2c_setscreen2 if expType D.
        %
        % so, expType C is always dSR3T2c_setscreens2
        if strcmpi(expType, 'C')
            trialTypes = repmat({'SR3'}, 1, nTrials);
        elseif strcmpi(expType, 'D')
            trialTypes = repmat({'SR4'}, 1, nTrials);
            if ptbParams.currSession <= 49
                % JerryLee sra3_2 before 090626 are actually SR3, but every sra3_2
                % after that is SR4. I haven't yet spotted a sure-fire way to know
                % that other than from the date.
                trialTypes = repmat({'SR3'}, 1, nTrials);
            elseif isfield(trials, 'SR3errorStrategy')
                has_err_strat = ~cellfun(@isempty, {trials.SR3errorStrategy});
                err_chng_to = [trials(has_err_strat).SR3errorStrategy];
                err_chng_to = {err_chng_to.error};
                % if we don't know the beginning error strategy, assume
                % resample.
                if ~has_err_strat(1)
                    has_err_strat(1) = true;
                    err_chng_to = ['resample' err_chng_to];
                end
                is_prob = false(nTrials, 1);
                err_chng_id = find(has_err_strat);
                for e_ix = 1:length(err_chng_id);
                    if e_ix == length(err_chng_id)
                        this_ix = err_chng_id(e_ix):nTrials;
                    else
                        this_ix = err_chng_id(e_ix):err_chng_id(e_ix+1)-1;
                    end
                    is_prob(this_ix) = strcmpi(err_chng_to{e_ix}, 'probation');
                end
                trialTypes(is_prob) = repmat({'SR3'}, 1, sum(is_prob));
            end
        end
    end
    
    %% Training levels for each trial.
    % Determines whether target,distractor were used.
    tl_colix = 1;
    if strcmpi(expType, 'D')
        tl_colix = 2;
    end
    trainingLevels = zeros(nTrials, 1);
    has_tl = false(nTrials, 1);
    %TODO: Should SR3trainingLevel or SR3InitialtrainingLevel get priority here?
    if isfield(trials, 'SR3trainingLevel')
        has_tl(~cellfun(@isempty, {trials.SR3trainingLevel})) = true;
        temp = cat(1, trials(has_tl & trialBool).SR3trainingLevel);
        trainingLevels(has_tl & trialBool) = temp(:, tl_colix);
    end
    if any(~has_tl & trialBool) && isfield(trials, 'SR3InitialtrainingLevel')
        temp = cat(1, trials(~has_tl & trialBool).SR3InitialtrainingLevel);
        trainingLevels(~has_tl & trialBool) = temp(:, tl_colix);
    end
    targBool = trainingLevels ~= 1;  % trial has target
    distBool = ~ismember(trainingLevels, [1 3 4]);  % trial has distractor
    clear tl_colix has_tl temp
    
    %     [trials(trialBool(~distBool)).newType] = deal('CentreOut');
    
    %% Get per-trial indices into targetTheta (targ/dist), cueColourMap, and per trial radii
    % method depends on trialType
    
    nExpTrials = sum(trialBool);
    targ_ix = nan(nExpTrials, 1);
    dist_ix = nan(nExpTrials, 1);
    cueCol_ix = nan(nExpTrials, 1);
    annulus_ix = nan(nExpTrials, 1);
    contexts = cell(nExpTrials, 1);
    
    for trialType = unique(trialTypes)
        typeBool = strcmpi(trialTypes, trialType);
        temp = classes(typeBool);
        switch trialType{1}
            case 'M'
                ptbParams.flipNames = {'fixationOnset' 'targetOnset' 'fixationOffset' 'imperativeCue' 'saccadeEnd' 'targetAcqLost'};
                %TODO: Check flipNames against the PTB code.
                %It was tough to follow, I don't know if the above is correct.
                % For "M" trials, we are looking for correct saccades typically after
                % flipScreen 2 or 3.
                % ??? Class 1 is up, then around the face of a clock in 45
                % degree increments until 8. If 9:16 present, those are the same except
                % larger amplitude.
                targ_ix(typeBool) = mod(temp, 8); targ_ix(targ_ix == 0) = 8;
                annulus_ix(typeBool) = 1 + (temp - targ_ix(typeBool))/8;
                cueCol_ix(typeBool) = ones(size(cueCol_ix(typeBool)));
            case 'A'
                warning('Not yet implemented expType A.');
                %trial class 0 or 1 only?
                
            case 'SR3'
                ptbParams.flipNames = {'fixationOnset' 'targetOnset' 'cueOnset' 'cueOffset' 'fixationOffset' 'targetAcqLost'};
                %imperativeCue == fixationOffset
                % trial class = 8*(cueRow-1) + 4*(cueIndex-1) + targetConfig
                % I rename targetConfig -> theta_targ_ix
                % cueRow is in 1:length(cueMatrix), chosen rand at the start of each block
                % targetConfig is in 1:4, representing 4 direction pairs, rand per block
                % cueIndex is in 1:2, when 1 the targ is in first half, 2: targ in
                % second half. targets are in [270 315 0 45 90 135 180 225]
                % targetChoice, in 1:8, = 4*(cueIndex-1) + targetConfig
                cueRow = ceil(temp/8);
                
                temp = temp - 8*(cueRow-1);
                cueIndex = ceil(temp/4);
                
                temp = temp - 4*(cueIndex-1);
                targetConfig = temp;
                
                targ_ix(typeBool) = 4*(cueIndex-1) + targetConfig;
                dist_ix(typeBool) = targ_ix(typeBool)...
                    + 4*double(targ_ix(typeBool)< 5)...
                    - 4*double(targ_ix(typeBool) > 4);
                
                cueCol_ix(typeBool) = sub2ind(size(cueMatrix), cueRow, cueIndex);
                annulus_ix(typeBool) = ones(size(annulus_ix(typeBool)));
                contexts(typeBool) = targetStr(targ_ix(typeBool));
            case 'SR4'
                ptbParams.flipNames = {'fixationOnset' 'targetOnset' 'cueOnset' 'cueOffset' 'fixationOffset' 'saccadeEnd' 'targetAcqLost'};
                % trial class = 108*(r-1) + 54*(y-1) + 9*(x-1) + z
                
                % r is cueColumn; 1 for u|l correct, 2 for d|r corr, rand per trial
                r = ceil(temp/108);
                temp = temp - 108*(r-1);
                
                % y is 1 for u/d, 2 for l/r, chosen randomly per block
                y = ceil(temp/54);
                temp = temp - 54*(y-1);
                
                % x is cueRow in cueMatrix, in 1:6, (determines targ/dist colour pairings per block)
                x = ceil(temp/9);
                temp = temp - 9*(x-1);
                
                % z in 1:9, rand per block. Indexes into target,distractor pairings
                % t order: 1 2 3 1 2 3 1 2 3 into targDistGroups(t_ix, r, y, 1)
                % d order: 1 1 1 2 2 2 3 3 3 into targDistGroups(d_ix, r, y, 2)
                z = temp;
                t_ix = mod(z, 3); t_ix(t_ix==0) = 3;
                d_ix = ceil(z/3);
                
                targ_ix(typeBool) = targDistGroups(sub2ind(size(targDistGroups), t_ix, r, y, 1*ones(size(t_ix))));
                dist_ix(typeBool) = targDistGroups(sub2ind(size(targDistGroups), d_ix, r, y, 2*ones(size(t_ix))));
                cueCol_ix(typeBool) = sub2ind(size(cueMatrix), x, r);
                annulus_ix(typeBool) = ones(size(annulus_ix(typeBool)));
                
                contexts(typeBool) = context_options((y-1)*2+r);
            otherwise
                error('trial type not recognized');
        end
    end
    
    
    %% Save to trial structure
    %.newType = 'M', 'A', 'SR3', or 'SR4'
    %.cueColour = 'r', 'g', or 'b'
    %.trainingLevel = 0: full exp; 1: fix only; 2: fix & dist; 3: full exp &
    %invis dist; 4: 3 + extra cue; 5: 4 + extra cue
    %.targRule = 'UU', 'UR', 'RR', etc. OR 'AnyUp', 'AnyRight', etc.
    %.targClass = in 1:8 for 8 locations, or in 1:16 with two annuli.
    %.targPol = [theta radius]
    %.targXY = [x y] coordinates in screen pixels.
    %(also distClass, distPol, distXY)
    
    [trials(trialBool).newType] = trialTypes{:};
    
    %Cue colours
    trCueColour = cueMatrix(cueCol_ix);
    trCueColour = cueColourMap(trCueColour);
    [trials(trialBool).cueColour] = trCueColour{:};
    
    %Training levels
    trainingLevels = num2cell(trainingLevels);
    [trials(trialBool).trainingLevel] = trainingLevels{:};
    
    %Contexts
    [trials(trialBool).targRule] = contexts{:};
    
    % Save target/distractor to trial structure
    trRadii = radius(annulus_ix);
    
    targDist = {'targ' 'dist'};
    td_ix = [targ_ix dist_ix];
    td_bool = [targBool distBool];
    for td = 1:2 %For targets and distractors.
        this_ix = td_ix(:, td);  % Indices into targetTheta
        this_bool = td_bool(:, td);  % If this trial had a targ/dist
        % stimulus class of 8 (or 16) possible locations
        trClass = nan(nTrials, 1);
        trClass(this_bool) = this_ix(this_bool);
        trClass(annulus_ix > 1) = annulus_ix(annulus_ix > 1).*trClass(annulus_ix > 1);
        trClass = num2cell(trClass);
        [trials(trialBool).([targDist{td} 'Class'])] = trClass{:};
        % Angle - needed for below
        trTheta = nan(nTrials, 1);
        trTheta(this_bool) = targetTheta(this_ix(this_bool));
        % Pol, assuming centre = 0, 0
        trPol = num2cell([trTheta trRadii], 2);
        [trials(trialBool).([targDist{td} 'Pol'])] = trPol{:};
        % Coordinates in screen pixels
        [trX, trY] = pol2cart(trTheta, trRadii);
        trXY = round(bsxfun(@plus, centrePt + fixJitterXY, [trX trY]));
        trXY = num2cell(trXY, 2);
        [trials(trialBool).([targDist{td} 'XY'])] = trXY{:};
        %String representing direction.
        trStr = cell(nTrials, 1);
        trStr(this_bool) = targetStr(this_ix(this_bool));
        [trials(trialBool).([targDist{td} 'Str'])] = trStr{:};
    end
    
end

%% Block
% The 'block' is considered to change whenever the trial type changes, the
% rule-pair changes, or the training-level changes long-enough to include a
% rule-pair change.

% We will store trials.newBlock and trials.hiatus
[trials.hiatus] = deal(false);

%Get the types, training levels, and rule ids.
[~, typ_id] = ismember({trials.newType}, unique({trials.newType}));
[~, trn_id] = ismember([trials.trainingLevel], unique([trials.trainingLevel]));
if isfield(trials, 'block')
    rul_id = [trials.block];
else
    rul_id = ones(1, nTrials);
    warning('TODO: Manually determine when the rule changes for sessions missing ptb.eData.trial.block');
end

%The trials when the type or rule changes are always block changes:
block_chng = diff([0 typ_id]) ~= 0 ...
           | diff([rul_id(1) rul_id]) ~= 0;

%When the training level changes may or may not be a block change
trlvl_chng = diff([trn_id(1) trn_id]) ~= 0;

%For each each trlvl_chng, determine if it was also a block change.
chng_ids = find(trlvl_chng);  %suspect 594 and 1262
for chng_ix = 1:length(chng_ids)
    if chng_ix < length(chng_ids)
        %If not the last, check to see if this is ending a hiatus
        is_end = trn_id(chng_ids(chng_ix)) == ...
            trn_id(find(block_chng(1:chng_ids(chng_ix)), 1, 'last'));
        if ~is_end
            %If it is not ending a hiatus, check to see if the block changes before next
            is_change = any(block_chng(chng_ids(chng_ix):chng_ids(chng_ix+1)-1));
            %if not ending hiatus, and block doesn't change, mark from here
            %'til next as part of the hiatus
            if ~is_change
                [trials(chng_ids(chng_ix):chng_ids(chng_ix+1)).hiatus] = deal(true);
            end
        end
    else
        %If the last training level change
        % 1 - Check to see if the block changes before the end
        is_change = any(block_chng(chng_ids(chng_ix):end));
        % 2 - If no more changes, check if this is not a return to
        % the training level used when the block started
        if ~is_change
            is_change = trn_id(chng_ids(chng_ix)) ~= ...
                trn_id(find(block_chng(1:chng_ids(chng_ix)), 1, 'last'));
        end
    end
    block_chng(chng_ids(chng_ix)) = is_change;
end

newBlock = num2cell(cumsum(block_chng));
[trials.newBlock] = newBlock{:};