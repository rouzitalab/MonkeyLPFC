function [IsInCondition] = getCornerConditionMembership(ptbTrials)
%[IsInCondition] = getCornerConditionMembership(ptbTrials)
%Determines whether or not each trial belongs to a certain condition.
%Considers only trials where the target was in a corner.
%IsInCondition is of size:
%nTrials x 4 corners x 2 (2 rules per corner) x 3 distractors x 3 colours
%
% Target locations are:
% 8 1 2
% 7   3
% 6 5 4
%
% Thus, corners are 2 4 6 8
%
% Each corner has two rules. e.g., 2 has AnyUp and AnyRight
% For each corner * rule, there are 3 distractors, ordered as near-mid-opp
% And then everything has 'r' 'g' 'b' cueColour.
%
% Examples:
%
% To get all trials to a corner with distractors in the opposite corner:
% squeeze(any(IsInCondition(:, corner_ix, :, 3, :), 5))
%
% To get all trials to a corner, using all distractors:
% squeeze(any(any(IsInCondition(:, corner_ix, :, :, :), 5), 4));

corners = [2 4 6 8];
rules = {{'AnyRight' 'AnyUp'}, {'AnyRight' 'AnyDown'}, {'AnyLeft' 'AnyDown'}, {'AnyLeft' 'AnyUp'}};
distractors = {{[8 7 6], [4 5 6]}, {[6 7 8], [2 1 8]}, {[4 3 2], [8 1 2]}, {[2 3 4], [6 5 4]}};
colours = {'r' 'g' 'b'};

nTrials = length(ptbTrials);
IsInCondition = false(nTrials, 4, 2, 3, 3);

for corner_ix = 1:4
    cornerBool = [ptbTrials.targClass] == corners(corner_ix);
    for rule_ix = 1:2
        ruleBool = strcmpi({ptbTrials.targRule}, rules{corner_ix}{rule_ix});
        for dist_ix = 1:3
            distBool = [ptbTrials.distClass] == distractors{corner_ix}{rule_ix}(dist_ix);
            for colour_ix = 1:3
                colourBool = strcmpi({ptbTrials.cueColour}, colours{colour_ix});
                
                IsInCondition(:, corner_ix, rule_ix, dist_ix, colour_ix) = ...
                    cornerBool & ruleBool & distBool & colourBool;
            end
        end
    end
end

