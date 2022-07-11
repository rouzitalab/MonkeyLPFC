function [mu, sigma] = normmle(X,ind)
%NORMMLE  Calculate the (biased) ML estimates for a normal pdf.
%   [MU, SIGMA] = NORMMLE(X) returns the sample means of X in the vector
%   MU, and returns the (biased) estimate of the covariance
%   matrix in SIGMA.
%
%   [MU, SIGMA] = NORMMLE(X,IND) is as above if IND is equal to 0,
%   otherwise (e.g., if IND = 1) SIGMA is a vector of variances
%   (i.e., here, we make an independence assumption).
%
%   Note that X contains M rows of sample points on N columns of
%   variables.

% Scott Gaffney   26 April 1998
% Department of Information and Computer Science
% University of California, Irvine

if (exist('ind','var') ~=1 | ind==0)
    ind = 0;
else
    ind = 1;
end

mu = mean(X)';
if (ind==1)
    sigma = var(X);
    sigma = sigma(:);
else
    sigma = cov(X,1); 
end
