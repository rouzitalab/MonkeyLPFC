function X = regmat(x,order)
%REGMAT  Make a Vandermonde regression matrix
%   X = REGMAT(x,order)


% Scott Gaffney   1 June 2002
% Department of Information and Computer Science
% University of California, Irvine


PROGNAME = 'regmat';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

x = x(:);
n = length(x);
X = [ones(n,1) x];
X = cumprod([X (X(:,2)*ones(1,order-1))],2);
