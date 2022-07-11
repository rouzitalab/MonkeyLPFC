function M = kmeans_capi(trajs,K,Ops)
%KMEANS_CAPI  Curve API for K-means.
%
%   Model = KMEANS_CAPI(Trajs,K,[Options])
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find

% Scott Gaffney   9 October 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'kmeans_capi';
METHOD = PROGNAME;
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

%%% Handle Argument Processing
%%%
n = nargin;
if (n<2)
  error([PROGNAME, ': incorrect number of parameters provided.']);
end
%%
%%% End Argument Processing

% extract the data
[X,D] = trajs2seq(trajs,Ops.zero,Ops.MinLen,'matrix');

% run K-means
[M.Mu,M.C] = kmeans(X,K);

% convert the output to "curve form"
M.Alpha = hist(M.C,(1:K)');
M.Alpha  = M.Alpha./sum(M.Alpha);
% put Mu into P,K,D form
if (D>1)
  P = size(X,2)/D;
  M.Mu = reshape(M.Mu,[P D K]);
  M.Mu = permute(M.Mu,[1 3 2]);
end
M.Options = Ops;
M.zero = Ops.zero;
M.method = 'kmeans';

% calculate the number of independent parameters
[P,K,D] = size(M.Mu);
M.NumIndParams = K*P*D;   %  mu
M.TrainLhood=[];


