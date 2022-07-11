function [Lhood,other] = lrm_cd_ab_like(varargin)
%LRM_CD_AB_LIKE  Calculate log-likelihood with LRM_CD_AB model.
%
%   [Lhood,Other] = LRM_CD_AB_LIKE(M,Trajs,[Options])
%    - M       : trained model
%    - Trajs   : 'Trajs' structure; See also CCTOOLBOX
%    - Options : see MODEL_LIKE
%
%   [Lhood,Other] = LRM_CD_AB_LIKE(M,x,Y,Seq,[Options])
%    - M       : trained model
%    - X,Y,Seq : curves in Sequence format; See also CCTOOLBOX
%    - Options : see MODEL_LIKE
%
%   Other
%     .C   : classification labels

% Scott Gaffney   15 February 1999
% DataLab@UCI
% Department of Information and Computer Science
% University of California, Irvine, USA.

PROGNAME = 'lrm_cd_ab_like';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%%% Handle Argument Processing
%%%
args = varargin; clear varargin;
n = length(args);
trajs=[]; x=[]; Y=[]; Seq=[]; Ops=[];
%
% Check for calling convention
%
% LRM_CD_AB_LIKE(M,Trajs,[Options])
if (n<4)
  M = args{1};
  trajs = args{2};
  if (n>2)
    Ops = args{3};
  end
  
% LRM_CD_AB_LIKE(M,x,Y,Seq,[Options])
else
  M = args{1};
  x = args{2};
  Y = args{3};
  Seq = args{4};
  if (n>4)
    Ops = args{5};
  end
end
%%
%%% End Argument Processing

% Handle Ops
if (~isfield(Ops,'NumSamps') | isempty(Ops.NumSamps))
  Ops.NumSamps = 500;
end
if (~isfield(Ops,'Scale') | isempty(Ops.Scale))
  Ops.Scale = 1.25;
end
if (~isfield(Ops,'VarScale') | isempty(Ops.VarScale))
  Ops.VarScale = 1;
end
if (~isfield(Ops,'MaxAtts') | isempty(Ops.MaxAtts))
  Ops.MaxAtts = 3;
end

% build data
if (isempty(Y))
  [Y,x,Seq] = trajs2seq(trajs,M.zero,M.Options.MinLen);
end

% Model setup
M.Mu = permute(M.Mu,[1 3 2]);
N = size(Y,1);

% calculate the likelihood
Pik = CalcPik(M,x,Y,Seq,Ops);
Lhood = sum(log(sum(Pik,2)))./prod(size(Y));
[trash, other.C] = max(Pik,[],2);

  
  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcPik 
%%
function pik = CalcPik(M,x,Y,Seq,Ops)
% Numerical integration
NumSamps = Ops.NumSamps;
MaxTries = Ops.MaxAtts;
Scale = Ops.Scale;
VarScale = Ops.VarScale;
[N,D] = size(Y);
n = length(Seq)-1;
P = M.order+1;
K = M.K;
Pid = zeros(n,D);
Pik = zeros(n,K);
R = M.R;  S = M.S;  T = M.T;  U = M.U;
TotalSamps = 0;
tries = 1;

while (1)
  TotalSamps = TotalSamps + NumSamps;

  % calculate the density at sampled points
  for k=1:K
    a = randn(NumSamps,1).*sqrt(R(k)) + 1;  % sample from N(1,r)
    b = randn(NumSamps,1).*sqrt(S(k));   % sample from N(0,s)
    for j=1:NumSamps
      Xhat = regmat(a(j)*x-b(j),P-1);
      XMu = Xhat*M.Mu(:,:,k);
      for d=1:D
        sigma = M.Sigma(k,d);  t = T(k,d);  u = U(k,d);
        for i=1:n
          indx = Seq(i):Seq(i+1)-1;  ni = length(indx);
          iR = eye(ni)./sigma - XMu(indx,d)*XMu(indx,d)'/(sigma^2/u ...
            + sigma*(XMu(indx,d)'*XMu(indx,d)));  siR = sum(iR);
          iV = iR - sum(iR,2)*siR/(1/t + sum(siR));
          Pid(i,d) = mvnormpdf_inv(Y(indx,d)',Xhat(indx,:)*M.Mu(:,d,k),iV);
        end
      end
      Pik(:,k) = Pik(:,k) + prod(Pid,2);
    end
  end
  pik = (Pik./TotalSamps) .* (ones(n,1)*M.Alpha'); % we keep Pik for next try!
  if (all(sum(pik,2))), break; end

  % we have detected some zeros, try again?
  if (tries==MaxTries)
    fprintf('lrm_aa_ta_sh_like: Integration failed, using realmin*1e100 instead.\n');
    zero = find(sum(pik,2)==0);
    pik(zero,:) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
    break;
  else
    fprintf('lrm_aa_ta_sh_like: Zero membership detected, trying integration again: %d\n',tries);
    tries = tries+1;
    S = VarScale*S;  % biased, but gets over some tricky integrations
    R = VarScale*R;
  end
end
