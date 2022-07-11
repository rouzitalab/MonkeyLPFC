function [Lhood,other] = lrm_cd_like(varargin)
%LRM_CD_LIKE  Calculate log-likelihood with LRM_CD model.
%
%   [Lhood,Other] = LRM_CD_LIKE(M,Trajs,[Options])
%    - M       : trained model
%    - Trajs   : 'Trajs' structure; See also CCTOOLBOX
%    - Options : see MODEL_LIKE
%
%   [Lhood,Other] = LRM_CD_LIKE(M,X,Y,Seq,[Options])
%    - M       : trained model
%    - X,Y,Seq : curves in Sequence format; See also CCTOOLBOX
%    - Options : see MODEL_LIKE
%
%   Other
%     .C   : classification labels

% Scott Gaffney   5 October 2003
% DataLab@UCI
% Department of Information and Computer Science
% University of California, Irvine, USA.

PROGNAME = 'lrm_cd_like';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%%% Handle Argument Processing
%%%
args = varargin; clear varargin;
n = length(args);
trajs=[]; X=[]; Y=[]; Seq=[]; Ops=[];
%
% Check for calling convention
%
% LRM_CD_LIKE(M,Trajs,[Options])
if (n<4)
  M = args{1};
  trajs = args{2};
  if (n>2)
    Ops = args{3};
  end
  
% LRM_CD_LIKE(M,X,Y,Seq,[Options])
else
  M = args{1};
  X = args{2};
  Y = args{3};
  Seq = args{4};
  if (n>4)
    Ops = args{5};
  end
end
%%
%%% End Argument Processing

if (isempty(Y))
  [Y,X,Seq] = trajs2seq(trajs,M.zero,M.Options.MinLen);
end
if (size(X,2)~=size(M.Mu,1))
  X = regmat(X,M.order);  
end


[N,D] = size(Y);
NumPoints = N*D;
K = M.K;
n = length(Seq)-1;
M.Mu = permute(M.Mu,[1 3 2]);

Pikd = zeros(n,K,D);
for k=1:K
  for d=1:D
    Mu       = M.Mu(:,d,k);
    sigma    = M.Sigma(k,d);
    r        = M.R(k,d);
    s        = M.S(k,d);
    for i=1:n
      indx   = Seq(i):Seq(i+1)-1;
      I      = eye(length(indx));
      XMu    = X(indx,:)*Mu;
      iR     = I/sigma - XMu*XMu'/(sigma^2/r + sigma*(XMu'*XMu));
      siR = sum(iR);
      iV = iR - sum(iR,2)*siR/(1/s + sum(siR));
      Pikd(i,k,d) = mvnormpdf_inv(Y(indx,d)',XMu',iV);
    end
  end
end
Pikd(:,:,1) = prod(Pikd,3);
Pikd(:,:,1) = Pikd(:,:,1) .* (ones(n,1)*M.Alpha');
s = sum(Pikd(:,:,1),2);
if (~all(s))
  fprintf([PROGNAME, ': log(0) detected, using log(K*realmin*1e100).\n']);
  zero = find(s==0);
  Pikd(zero,:,1) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
  s(zero) = sum(Pikd(zero,:,1),2);
end
Lhood = sum(log(s));
Lhood = Lhood./NumPoints;

% Classify the sequences
if (nargout>1)
  [trash, other.C] = max(Pikd(:,:,1),[],2);    
end
