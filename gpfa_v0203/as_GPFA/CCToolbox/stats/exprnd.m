function A = exprnd(mu,m,n);
%EXPRND Sample from exponential distribution.
%   A = EXPRND(MU) is a random matrix from the exponential distribution 
%   with mean MU (A is always the same size as MU).
%
%   A = EXPRND(MU,M,N) is as above but A is M-by-N. 

% Scott J Gaffney   5 February 2003
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'exprnd';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

m = cexist('m',size(mu,1));
n = cexist('n',size(mu,2));

A = zeros(m, n);
A = -mu .* log(rand(m,n));

anyNaN = find(mu<=0);
if (~isempty(anyNaN))
  if (prod(size(mu))==1)
    A = NaN + ones(m,n);
  else
    A(anyNaN) = NaN;
  end
end

