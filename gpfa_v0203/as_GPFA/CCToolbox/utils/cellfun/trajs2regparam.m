function B = trajs2regparam(trajs,order,zero)
%trajs2regparam  Convert cell array of trajs to regression parameters
%   B = trajs2regparam(TRAJS,order,['zero']) fits regression lines to each
%   curve in TRAJS and returns them in B.
%
%   B = trajs2regparam(Y,order,['zero'])
%

% Scott J Gaffney   10 December 2002
% Department of Information and Computer Science
% University of California, Irvine.
%
% Changes
% ---------------------------------
%

PROGNAME = 'trajs2regparam';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

% args
DoZero = 0;
if (exist('zero')==1 & isstr(zero) & strcmp(zero,'zero'))
  DoZero = 1;
end


if (iscell(trajs))
  n = length(trajs);
  [mnlen, lens] = meanlength(trajs);
  maxlen = max(lens);
  D = size(trajs{1},2);
  x = regmat((0:maxlen-1)',order);
  B = zeros(n,D*(order+1));
  for j=1:n
    y = trajs{j};
    if (DoZero)
      y = y - ones(size(y,1),1)*y(1,:);
    end
    tmp = wls(x(1:size(y,1),:),y);
    B(j,:) = tmp(:)';
  end
  
else
  [n,m] = size(trajs);
  x = regmat((0:m-1)',order);
  B = zeros(n,order+1);
  for j=1:n
    y = trajs(j,:)';
    if (DoZero)
      y = y - y(1);
    end
    tmp = wls(x,y);
    B(j,:) = tmp(:)';
  end
end
B = B';