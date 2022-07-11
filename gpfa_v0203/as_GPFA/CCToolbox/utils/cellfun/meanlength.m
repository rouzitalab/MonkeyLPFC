function [mn,len] = meanlength(trajs)
%MEANLENGTH  Calculate mean length over all the cells.
%   [MEAN,LENGTHS] = MEANLENGTH returns the mean lenght over all cells
%   in MEAN and also a vector of lengths in LENGTHS.

% Scott J Gaffney   11 September 2001
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'meanlength';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


for i=1:length(trajs)
  len(i) = size(trajs{i},1);
end

len=len(:);
mn = mean(len);
