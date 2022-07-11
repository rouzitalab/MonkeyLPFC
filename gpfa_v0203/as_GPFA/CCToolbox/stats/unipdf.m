function y = unipdf(X,a,b)
%UNIPDF  Calculate univariate uniform density values at X.
%   UNIPDF(X,a,b) is a vector of density values calculated at X from
%   the uniform density over the interval [a,b].
%
%   Note that X contains N rows of sample points, 'a','b' are scalars.

% Scott Gaffney   21 November 2001
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'unipdf';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

a = cexist('a',0);
b = cexist('b',1);

X=X(:);
a = a(:);
b = b(:);

if (length(a)~=1 & length(a)~=length(X))
  error([PROGNAME, ': the length of a is invalid.']);
end
if (length(b)~=1 & length(b)~=length(X))
  error([PROGNAME, ': the length of b is invalid.']);
end

where = find(X>=a & X<=b);
y = zeros(length(X),1);
y(where) = 1./(b-a);