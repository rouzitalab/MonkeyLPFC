function [b,sigma] = wls(X,Y,W,NoNorm)
%WLS  Weighted Least Squares (WLS).
%   WLS(X,Y,W) is the vector of regression coefficients resulting
%   from the weighted regression of Y on X. W contains the weights.
%
%   [b,sigma_hat] = WLS(...,NoNorm) will return the estimated variance
%   of the error term if you include a second output parameter. 
%   However, sigma_hat is normalized by sum(W) if NoNorm does not
%   exist or if it is not equal to 1.
%
%   [...] = WLS(X,Y) assumes W==ones(size(X,1),1).

% Scott Gaffney   2 October 1998
% Department of Information and Computer Science
% University of California, Irvine

% References:
%  [1] N. Draper and H. Smith (1981). "Applied Regression Analysis, Second
%       Edition." New York: Wiley.

% Changes
% -------------------------------
% 24 September 2001 (Scott Gaffney)
%   Changed the memory hog by replacing the SPARSE(DIAG(W)) call
%   with SPDIAGS.

PROGNAME = 'wls';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


[n,p] = size(X);
[yn,D] = size(Y);
if (n~=yn)
    error('WLS: the number of rows in X must equal those in Y.');
end
if (exist('W','var') ~= 1)
  W = ones(n,1);
end
W = W(:);
if (min(size(W)) ~= 1)
    error ('WLS: W must be a vector.');
end

% Construct the inverse covariance matrix inv(V) which consists of
% the elements of W along the diagonal. See [1], pg. 111.
% Vinv = spdiags(W,0,length(W),length(W));  % sparse diagonal matrix
% rtVinv = sqrt(Vinv);  % we need the square root to premultiply

% Find the least squares solution.
% [QY,R] = qr(rtVinv*X,rtVinv*Y);  % for sparse matrices only
% b = R\QY;

% [Q, R] = qr(rtVinv*X,0);  % Use this instead of qr(rtVinv*X)!!
% b = R\(Q'* rtVinv*Y);

if (D>p), rtWmat = sqrt(W)*ones(1,D);  else rtWmat = sqrt(W)*ones(1,p);  end
VX = rtWmat(:,1:p).*X;
VY = rtWmat(:,1:D).*Y;
[Q,R] = qr(VX,0);
b = R\(Q'*VY);

% Return the estimate for the weighted variance of the regression.
if (nargout > 1)
  if (exist('NoNorm') ~= 1 | NoNorm ~= 1)
    df = sum(W);     % ML
    % df = sum(W)-p;  % non-biased
  else
    df = 1;
  end
  residuals = (Y-X*b);
  Wmat = W*ones(1,D);
  sigma =  ((Wmat.*residuals)'*residuals)./df;
%   sigma =  (residuals'*Vinv*residuals)./df;
end
clear Wmat;
