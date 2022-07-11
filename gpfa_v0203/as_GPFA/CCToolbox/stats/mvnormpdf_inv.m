function y = mvnormpdf_inv(X,mu,isigma)
%MVNORMPDF_INV  Calculate multivariate normal density values at X.
%   MVNORMPDF_INV(X,MU,iSigma) is a vector of density values calculated at X from
%   the multivariate normal density with parameters MU, inv(iSigma). You must
%   pass-in the inverse of the covariance matrix to this function. See
%   MVNORMPDF if you want to pass-in the simple covariance matrix instead.
%
%   Note that X contains N rows of sample points, each of dimension
%   D (i.e., size(X)==[N,D]), MU is either a N by D matrix giving
%   N different means for each of the random variables in X, or a
%   vector of length D functioning as the common mean for all N random
%   variables, and iSigma is always a full D by D inverse covariance matrix.
%   See MVNORMPDF if you have a special covariance structure.

% Scott Gaffney   10 June 1998
% DataLab @UCI
% Department of Information and Computer Science
% University of California, Irvine
%
% SJG: 14 February 1999
% added the ability to provide a vector sized sigma argument.
%
% 25 September 2001 (Scott Gaffney)
%   Replaced DIAG with SPDIAGS to avoid possible memory blow-ups.
%
% 22 May 2003 (Scott Gaffney)
%   Specialized MVNORMPDF to this file MVNORMPDF_INV.

PROGNAME = 'mvnormpdf_inv';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

sz_mu  = size(mu);
psz_mu = prod(sz_mu);
[n,d]  = size(X);

% use common mean?
if (psz_mu~=n*d)
  mu = ones(n,1)*mu(:)';  % ...then copy the mean

% If X is a row vector, then make sure mu is also
elseif (n==1)     % isRowVector?
  mu = mu(:)';    % make row vector
end


dsigma = 1/det(isigma);

y = zeros(n,1);    % density outputs
X = X-mu;          % centered data
y = exp(-0.5 .* sum(X*isigma.*X,2)) ./ ((2*pi).^d *dsigma).^(.5);

