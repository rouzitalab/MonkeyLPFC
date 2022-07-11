function y = mvnormpdf(X,mu,sigma)
%MVNORMPDF  Calculate multivariate normal density values at X.
%   MVNORMPDF(X,MU,SIGMA) is a vector of density values calculated at X from
%   the multivariate normal density with parameters MU, SIGMA.
%
%   Note that X contains N rows of sample points, each of dimension
%   D (i.e., size(X)==[N,D]), MU is either a N by D matrix giving
%   N different means for each of the random variables in X, or a
%   vector of length D functioning as the common mean for all N random
%   variables, and SIGMA is either a single variance value, a D by 1 vector 
%   of variances, or a full D by D covariance matrix.

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
%   optimized preprocessing

PROGNAME = 'mvnormpdf';
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
elseif (n==1)  % isRowVector?
  mu = mu(:)';       % make row vector
end


sz = size(sigma);
FULL_LIM = 1000;
% Deal with scalar and vector separately
if (min(sz)==1) % isScalarOrVector?

  % scalar
  if (max(sz)==1)
    dsigma = sigma.^d;
    sigma = 1/sigma;
    if (d > FULL_LIM)
      isigma = spdiags(sigma(ones(d,1)),0,d,d);
    else
      isigma = diag(sigma(ones(d,1)));
    end
    
  % vector
  else
    dsigma = prod(sigma);
    if (d > FULL_LIM)
      isigma = spdiags(1./sigma,0,d,d);
    else
      isigma = diag(1./sigma);
    end
  end
  
% Full covariance
else
  isigma = inv(sigma);
  dsigma = det(sigma);
end
clear sigma;

%y = zeros(n,1);  % density outputs
X = X-mu;                % centered data
y = exp(-0.5 .* sum(X*isigma.*X,2)) ./ ((2*pi).^d * dsigma).^(.5);
if (d > FULL_LIM)
  y = full(y);
end

