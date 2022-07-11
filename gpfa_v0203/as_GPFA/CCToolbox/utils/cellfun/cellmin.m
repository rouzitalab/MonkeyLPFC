function mn = cellmin(x)
%CELLMIN  Calculate min of each cell
%   MIN = CELLMIN(X)

% Scott J Gaffney   11 September 2001
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'cellmin';
if (~nargin)
  help(PROGNAME);
  return;
end


sz = size(x);
if (sz(1)~=1 & sz(2)~=1)
  error('CELLMIN only supports vector-sized cell arrays');
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

for i=1:length(x)
  if (~isempty(x{i}))
    mn(i,:) = min(x{i},[],1);
  end
end

