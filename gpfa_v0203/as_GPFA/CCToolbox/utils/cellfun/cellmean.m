function mn = cellmean(x)
%CELLMEAN  Calculate mean of each cell
%   X = CELLMEAN(X)

% Scott J Gaffney   11 September 2001
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'cellmean';
if (~nargin)
  help(PROGNAME);
  return;
end

sz = size(x);
if (sz(1)~=1 & sz(2)~=1)
  error('CELLMEAN only supports vector-sized cell arrays');
end

% set the size of the output array
mn=[];
a = nan;
for i=1:length(x)
  if (~isempty(x{i}))
    mn = a(ones(length(x),size(x{i},2))); % set to NaNs
    break;
  end
end
if (isempty(mn))  % if empty then just return
  return;
end

% otherwise, calculate the means
for i=1:length(x)
  if (~isempty(x{i}))
    mn(i,:) = mean(x{i},1);
  end
end

