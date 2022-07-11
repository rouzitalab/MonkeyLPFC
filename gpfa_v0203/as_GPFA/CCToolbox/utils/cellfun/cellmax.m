function mx = cellmax(x)
%CELLMAX  Calculate max of each cell
%   MAX = CELLMAX(X)

% Scott J Gaffney   11 September 2001
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'cellmax';
if (~nargin)
  help(PROGNAME);
  return;
end


sz = size(x);
if (sz(1)~=1 & sz(2)~=1)
  error('CELLMAX only supports vector-sized cell arrays');
end

% set the size of the output array
mx=[];
a = nan;
for i=1:length(x)
  if (~isempty(x{i}))
    mx = a(ones(length(x),size(x{i},2))); % set to NaNs
    break;
  end
end
if (isempty(mx))  % if empty then just return
  return;
end

for i=1:length(x)
  if (~isempty(x{i}))
    mx(i,:) = max(x{i},[],1);
  end
end

