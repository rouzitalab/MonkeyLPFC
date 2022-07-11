function C = kmeans_class(M,trajs)
%KMEANS_CLASS  Classify test data
%   C = KMEANS_CLASS(Model,trajs)

% Scott J. Gaffney	10 November 2004
% Department of Information and Computer Science
% University of California, Irvine


PROGNAME = 'kmeans_class';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


[ni,K,D] = size(M.Mu);
Y = trajs2seq(trajs,M.zero,M.Options.MinLen,'matrix');
[n,P] = size(Y);

mu_PxK = reshape(permute(M.Mu,[1 3 2]),[ni*D K]);
mu_nxPxK = permute(mu_PxK(:,:,ones(n,1)),[3 1 2]);
Y_nxPxK = Y(:,:,ones(1,K));

Eik = permute(sum((Y_nxPxK - mu_nxPxK).^2,2),[1 3 2]);
[trash C] = min(Eik,[],2);  % assign new labels
  
