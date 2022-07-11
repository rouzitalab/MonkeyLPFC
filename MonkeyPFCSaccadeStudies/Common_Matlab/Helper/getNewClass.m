function [trials] = getNewClass(trials, varargin)
%[trials] = getNewClass(trials)
%adds a field .newClass to trials structure.
%
%[trials] = getNewClass(trials, params)
% params is a structure with field classifyTarg
%For 'class', 'sacClass', 'targClass', .newClass simply takes that value.
%For 'cueColour', .newClass is the string colour of the cue (r,g,b)
%For 'sacPol', 'targPol', .newClass is the direction string (UU,UR,RR... etc)
%For 'regionRule', .newClass is the cued rule ('upper', 'lower', 'rightward',
%'leftward'); only works for the RegionRule experiment (SR3_2)

default_params.classifyTarg = 'class';
params = varg2params(varargin, default_params);

if ismember(params.classifyTarg, {'class', 'sacClass', 'targClass'})
    newClass = [trials.(params.classifyTarg)];
    newClassStr = arrayfun(@num2str, newClass, 'unif', 0);
elseif strcmpi(params.classifyTarg, 'cueColour')
    classStr = {'r' 'g' 'b'};
    [~, newClass] = ismember({trials.cueColour}, classStr);
    newClassStr = classStr(newClass);
elseif ismember(params.classifyTarg, {'sacPol', 'targPol'})
    
%     [~,ia,~] = unique([trials.targClass]);
%     targetTheta = cat(1, trials(ia).targPol);
%     targetTheta = targetTheta(:, 1);
    
    targetTheta = [3*pi/2 7*pi/4 0 pi/4 pi/2 3*pi/4 pi 5*pi/4];
    targetStr = {'UU' 'UR' 'RR' 'DR' 'DD' 'DL' 'LL' 'UL'};
    
    %histc requires increasing values, need to sort, but remember order
    [tmpTargetTheta, ix] = sort(targetTheta);
    targetDiff = min(diff(tmpTargetTheta));
    tmpEdgesTheta = [tmpTargetTheta - targetDiff/2, tmpTargetTheta(end) + targetDiff/2];
    
    %Get the angles
    trPol = cat(1, trials.(params.classifyTarg));
    trPol = trPol(:, 1);
    trPol(trPol < min(tmpEdgesTheta)) = trPol(trPol < min(tmpEdgesTheta)) + 2*pi;
    
    %Bin the angles.
    [~, newClass] = histc(trPol, tmpEdgesTheta);
    newClass(isnan(trPol)) = nan;
    newClass(~isnan(newClass)) = ix(newClass(~isnan(newClass))); % undo sorting
    
    newClassStr = cell(size(newClass));
    newClassStr(~isnan(newClass)) = targetStr(newClass(~isnan(newClass)));
    
elseif strcmpi(params.classifyTarg, 'regionRule')
    classStr = {'AnyUp' 'AnyDown' 'AnyLeft' 'AnyRight'};
    [~, newClass] = ismember({trials.targRule}, classStr);
    newClassStr = classStr(newClass);
    
elseif strcmpi(params.classifyTarg, 'combinedRegionRule')
    %36 classes: 4 regions * 3 targs per region * 3 colours
    regStr = {'AnyUp' 'AnyDown' 'AnyLeft' 'AnyRight'};
    targetStr = {'UU' 'UR' 'RR' 'DR' 'DD' 'DL' 'LL' 'UL'};
    colStr = {'r' 'g' 'b'};
    reg_targ_map = [... %clock-wise within each region
        8 1 2;... 'AnyUp': UL, UU, UR
        4 5 6;... 'AnyDown': DR, DD, DL
        6 7 8;... 'AnyLeft': DL, LL, UL
        2 3 4];%  'AnyRight': UR, RR, DR
    
    [~, regClass] = ismember({trials.targRule}, regStr);
    trRegStr = regStr(regClass);
    [~, targetClass] = ismember({trials.targStr}, targetStr);
    trTargetStr = targetStr(targetClass);
    [~, colClass] = ismember({trials.cueColour}, colStr);
    trColStr = colStr(colClass);
    
    newClassStr = cell(1, length(trials));
    newClass = nan(1, length(trials));
    for tr_ix = 1:length(trials)
        newClassStr{tr_ix} = [trRegStr{tr_ix} '_' trTargetStr{tr_ix} '_' trColStr{tr_ix}];
        newClass(tr_ix) = find(reg_targ_map(regClass(tr_ix), :) == targetClass(tr_ix));
    end
    newClass = 9*(regClass-1) + 3*(newClass-1) + (colClass-1); %Base-zero
end

newClass = num2cell(newClass);
[trials.newClass] = newClass{:};
[trials.newClassStr] = newClassStr{:};