function v = cell2vector(C)
%CELL2VECTOR  Make a cell array into one long column vector
%   V = CELL2VECTOR(C)

% Scott J Gaffney   31 January 2002
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'cell2vector';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


if (~iscell(C))
  v = C;
  return;
end

v = [];
C = C(:);
len = length(C);
for (i=1:len)
  num = prod(size(C{i}));
  v(end+1:end+num,1) = C{i}(:);
end