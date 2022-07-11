function [LHOOD,C] = gmix_like(M,trajs)
%GMIX_LIKE  Calculate log-likelihood for Gaussian mixtures.
%
%   [Lhood,C] = GMIX_LIKE(M,X)
%    - M       : trained model
%    - X       : N rows of feature vectors
%
%   [Lhood,C] = GMIX_LIKE(M,Trajs])
%    - M       : trained model
%    - Trajs   : 'Trajs' structure; See also CCTOOLBOX
%    - Options : .zero

% Scott Gaffney   15 February 1999
% DataLab@UCI
% Department of Information and Computer Science
% University of California, Irvine, USA.

PROGNAME = 'gmix_like';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

% Handle preprocessing 
% if (length(size(M.Sigma))==2)
%   if (isvector(M.Sigma))
%     covType = 3;  % uniform variances
%   else
%     covType = 2;  % zero covariances
%   end
% else
%   covType = 1;  % full covariance matrices
% end
covType = 1;

% Handle invalid data
if (isnan(M.TrainLhood))
  LHOOD = NaN;
  C = ones(length(trajs),1);
  return;
end

% convert the data
X = trajs2seq(trajs,M.zero,M.Options.MinLen,'matrix');
  
[P,K,D] = size(M.Mu);
Pikd = zeros(size(X,1),K,D);
% Calculate the probability of the data
if (covType==3)
  for k = 1:K
    for d=1:D
      start = (d-1)*P+1;  stop = start+P;
      indx = start:stop-1;
      Pikd(:,k,d) = mvnormpdf(X(:,indx),M.Mu(:,k),M.Sigma(k));
    end
  end
elseif (covType==2)
  for k = 1:K
    for d=1:D
      start = (d-1)*P+1;  stop = start+P;
      indx = start:stop-1;
      Pikd(:,k,d) = mvnormpdf(X(:,indx),M.Mu(:,k),M.Sigma(:,k));
    end
  end
else
  for k = 1:K
    for d=1:D
      indx = (d-1)*P+1:d*P;  % un-concatenate the dimensions
      Pikd(:,k,d) = mvnormpdf(X(:,indx),M.Mu(:,k,d),M.Sigma(:,:,k,d));
    end
  end
end
Pik = prod(Pikd,3).*(ones(size(X,1),1)*M.Alpha');
LHOOD = sum(log(sum(Pik,2)))./prod(size(X));

if (nargout>1)
  [trash,C] = max(Pik,[],2);
  C = C(:);
end
