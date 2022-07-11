function y = normpdf(X,mu,Var)
%NORMPDF  Calculate univariate normal density values at X.
%   NORMPDF(X,MU,VAR) is a vector of density values calculated at X from
%   the normal density with parameters MU, VAR.
%
%   Note that X contains N rows of sample points, MU is either an
%   N by 1 vector giving N different means for each of the random 
%   variables in X, or a scalar functioning as the common mean for 
%   all N random variables. VAR is either a single variance value, or 
%   an N by 1 vector of variances for each random variable.

% Scott Gaffney   21 November 2001
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'normpdf';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

mu = cexist('mu',0);
Var = cexist('Var',1);

X=X(:);
mu = mu(:);
Var = Var(:);

if (length(mu)~=1 & length(mu)~=length(X))
  error([PROGNAME, ': the length of MU is invalid.']);
end
if (length(Var)~=1 & length(Var)~=length(X))
  error([PROGNAME, ': the length of VAR is invalid.']);
end

y = 1./(sqrt(2*pi*Var)) .* exp((-1/2)*(((X-mu).^2)./Var));




