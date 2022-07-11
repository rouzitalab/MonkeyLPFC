function [Lhood,other] = srm_ab_like(varargin)
%SRM_AB_LIKE  Calculate log-likelihood with SRM_AB model.
%
%   [Lhood,Other] = SRM_AB_LIKE(M,Trajs,[Options])
%    - M       : trained model
%    - Trajs   : 'Trajs' structure; See also CCTOOLBOX
%    - Options : see MODEL_LIKE
%
%   [Lhood,Other] = SRM_AB_LIKE(M,X,Y,Seq,[Options])
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

PROGNAME = 'srm_ab_like';
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
% SRM_AB_LIKE(M,Trajs,[Options])
if (n<4)
  M = args{1};
  trajs = args{2};
  if (n>2)
    Ops = args{3};
  end
  
% SRM_AB_LIKE(M,x,Y,Seq,[Options])
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

% Ops.MsgHnd = msgbar([],'');
% calculate the likelihood
[Pik, scale] = CalcPik(M,x,Y,Seq,Ops);
Lhood = (sum(log(sum(Pik,2))) + N*log(scale))./prod(size(Y));
[trash, other.C] = max(Pik,[],2);
% delete(Ops.MsgHnd);
  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcPik 
%%
function [Pik,scale] = CalcPik(M,x,Y,Seq,Ops)
% Numerical integration
NumSamps = Ops.NumSamps;
MaxTries = Ops.MaxAtts;
Scale = Ops.Scale;
VarScale = Ops.VarScale;

[N,D] = size(Y);
n = length(Seq)-1;
K = M.K;
mlen = max(diff(Seq));
Piid = zeros(N,D);
Piimk = zeros(N,NumSamps,K);
Pik = zeros(n,K);
S = M.S;  R = M.R;
TotalSamps = 0;
tries = 1;

while (1)
  TotalSamps = TotalSamps + NumSamps;
  Pik(:) = 0;

  % calculate the density at sampled points
  for k=1:K
    r = R(k);   s = S(k);
    a = randn(NumSamps,1).*sqrt(r) + 1;  % sample from N(1,r)
    b = randn(NumSamps,1).*sqrt(s);      % sample from N(0,s)
    for j=1:NumSamps
      Xhat = bsplinebasis(M.knots,M.order,a(j)*x-b(j));
      %if (mod(j,50)==0) msgbar(Ops.MsgHnd,sprintf('K %d, Sample %d\n',k,j)); end
      for d=1:D
        Piid(:,d) = normpdf(Y(:,d),Xhat*M.Mu(:,d,k),M.Sigma(k,d));
      end
      Piimk(:,j,k) = prod(Piid,2);
    end
  end
  
  % now scale the data to avoid underflow with long curves
  % and sum across the sample integration points
  scale = mean(mean(mean(Piimk)));
  Piimk_scl = Piimk./scale;
  for k=1:K
    for j=1:TotalSamps
      Pik(:,k) = Pik(:,k) + sprod(Piimk_scl(:,j,k),Seq,mlen);
    end
  end
  clear Piimk_scl;
  Pik = (Pik./TotalSamps) .* (ones(n,1)*M.Alpha');
  if (all(sum(Pik,2))), break; end

  % we have detected some zeros, try again?
  if (tries==MaxTries)
    fprintf(['srm_tt_sh_like: Integration failed, using realmin*1e100 ',...
      'instead.\n']);
    zero = find(sum(Pik,2)==0);
    Pik(zero,:) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
    break;
  else
    fprintf(['srm_tt_sh_like: Zero membership detected, trying ', ...
        'integration again: %d\n'],tries);
    tries = tries+1;
    S = VarScale*S;  % biased, but gets over some tricky integrations
    R = VarScale*R;  % biased, but gets over some tricky integrations
    Piimk = [zeros(N,NumSamps,K) Piimk]; % save current values
  end
end
