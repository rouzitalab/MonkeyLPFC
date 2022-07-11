function M = gmix_capi(trajs,K,Ops)
%GMIX_CAPI  Curve API for Gaussian mixtures.
%
%   Model = GMIX_CAPI(Trajs,K,[Options])
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find

% Scott Gaffney   9 October 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'gmix_capi';
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
Ops = cexist('Ops',[]);
%%
%%% End Argument Processing


% Preprocessing
Ops = DefaultOptions(Ops);

% extract the data
[X,D] = trajs2seq(trajs,Ops.zero,Ops.MinLen,'matrix');

% run Gaussian mixtures
M = gmix(X,K,Ops);

% convert the output to "curve form"
if (D>1)
  m = size(X,2);
  P = m/D;
  M.Mu = reshape(M.Mu,[P D K]);
  M.Mu = permute(M.Mu,[1 3 2]);
  M.concatSigma = M.Sigma;  M.Sigma=[];
  for (d=1:D)
    indx = (d-1)*P+1:d*P;
    M.Sigma(:,:,:,d) = M.concatSigma(indx,indx,:);
  end
  if (M.Options.Sigma.Share==1), shK = 1; else  shK=K; end
  if (M.Options.Sigma.Diagonal==0)
    M.NumIndParams = M.NumIndParams - shK*m*(m-P)/2;
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%

function Ops = DefaultOptions(Ops);
Ops = SetFieldDef(Ops,'zero','nozero');
Ops = SetFieldDef(Ops,'MinLen',[]);
