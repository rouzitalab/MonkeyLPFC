function x = normalize(x,dim)
%NORMALIZE  Normalize vector or matrix entries.
%   X = Normalize(X,dim)

% Scott Gaffney   11 January 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'normalize';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

dim = cexist('dim',1);
if (ndims(x)>2)
  error([PROGNAME,': error, you can only normalize 2-D matrices.\n']);
end

if (isvector(x))
  x = x/sum(x);
elseif (dim==1)
  x = x./ (ones(size(x,dim),1)*sum(x,dim)); % much quicker than repmat
elseif (dim==2)
  x = x./ (sum(x,dim)*ones(1,size(x,dim))); % much quicker than repmat
else
  error([PROGNAME,': error, you can only normalize across one of the', ...
        ' first two dimensions.\n']);
end