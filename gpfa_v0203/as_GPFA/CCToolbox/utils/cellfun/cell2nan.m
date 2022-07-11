function m = cell2nan(C)
%CELL2NAN  Make a cell array into a matrix with NaNs for missing values
%   M = CELL2NAN(C) only works for cell arrays with matrix elements.
%   If mlen = max(cell_len(C)) and D = size(C{1},2), then M is
%   of size n-by-mlen-by-D where n is length(C(:)).

% Scott J Gaffney   31 January 2002
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'cell2nan';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


if (~iscell(C))
  m = C;
  return;
end

C = C(:);
n = length(C);
mlen = max(cell_len(C));
D = size(C{1},2);
m = ones(mlen,D,n)*nan;

for (i=1:n)
  [r,c] = size(C{i});
  m(1:r,1:c,i) = C{i};
end
m = permute(m,[3 1 2]);