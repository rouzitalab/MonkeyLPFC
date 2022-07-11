function bool = isvector(x)
%ISVECTOR  Determine if X is a vector.
%   ISVECTOR(X)

% Scott Gaffney   11 January 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'isvector';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

bool=0;
if (prod(size(x))==length(x))
  bool = 1;
end