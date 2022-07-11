function lens = cell_len(trajs)
%CELL_LEN
%   lens = CELL_LEN

% Scott J Gaffney   20 May 2004
% School of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'cell_len';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


for i=1:length(trajs)
  lens(i) = size(trajs{i},1);
end
lens=lens(:);
