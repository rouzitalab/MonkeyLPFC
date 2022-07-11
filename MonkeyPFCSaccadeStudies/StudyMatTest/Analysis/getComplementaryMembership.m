function [PairBool] = getComplementaryMembership(ptbTrials)
%function [PairBool] = getComplementaries(ptbTrials)
% Finds pairs of complementary trials
% PairBool is of size nTrials, 4 target-pairs, 2 indices in pair, 3 colours
% You'll want to find trials with the same target-pair and same colour, but
% different indices in the pair
% do squeeze(sum(PairBool))
myColours = {'r' 'g' 'b'};
nTrials = length(ptbTrials);
PairBool = false(nTrials, 4, 2, 3);
for targ_ix = 1:4
    for pair_ix = 1:2
        targBool = [ptbTrials.sacClass] == targ_ix + (pair_ix-1)*4;
        for colour_ix = 1:3
            colourBool = strcmpi({ptbTrials.cueColour}, myColours{colour_ix});
            
            PairBool(targBool & colourBool, targ_ix, pair_ix, colour_ix) = true;
        end
    end
end