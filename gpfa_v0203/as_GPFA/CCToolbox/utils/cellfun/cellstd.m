function st = cellstd(x)
%CELLSTD  Calculate standard deviation of each cell
%   X = CELLSTD(X)

% Scott J Gaffney   11 September 2001
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'cellstd';
if (~nargin)
  help(PROGNAME);
  return;
end

sz = size(x);
if (sz(1)~=1 & sz(2)~=1)
  error('CELLSTD only supports vector-sized cell arrays');
end

% set the size of the output array
st=[];
a = nan;
for i=1:length(x)
  if (~isempty(x{i}))
    st = a(ones(length(x),size(x{i},2))); % set to NaNs
    break;
  end
end
if (isempty(st))  % if empty then just return
  return;
end

% otherwise, calculate the means
for i=1:length(x)
  if (~isempty(x{i}))
    st(i,:) = std(x{i});
  end
end

