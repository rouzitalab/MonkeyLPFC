function [b,sigma] = ols(X,Y,Order)
%OLS  Ordinary Least Squares (OLS).
%   OLS(X,Y,[ORDER]) is the vector of regression coefficients resulting
%   from the  regression of Y on X. If X's first column does not
%   contain a vector of ones, then it is added. If ORDER is provided
%   then X will be set up to facilitate regression of order ORDER.
%
%   [b,sigma_hat] = OLS(...) will return the estimated variance
%   of the error term if you include a second output parameter. 

% Scott Gaffney   9 October 2001
% Department of Information and Computer Science
% University of California, Irvine


PROGNAME = 'ols';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


[n,p] = size(X);

if (any(ones(n,1) ~= X(:,1)))
  X = [ones(n,1) X];
end

if (exist('Order')==1 & ~isempty(Order) & Order>1)
  X = cumprod([X repmat(X(:,2),1,Order-1)],2);
end

[b,sigma] = wls(X,Y);