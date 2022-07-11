function [Lhood,other] = lrm_b_like(varargin)
%LRM_B_LIKE  Calculate log-likelihood with LRM_B model.
%
%   [Lhood,Other] = LRM_B_LIKE(M,Trajs,[Options])
%    - M       : trained model
%    - Trajs   : 'Trajs' structure; See also CCTOOLBOX
%    - Options : see MODEL_LIKE
%
%   [Lhood,Other] = LRM_B_LIKE(M,x,Y,Seq,[Options])
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

PROGNAME = 'lrm_b_like';
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
% LRM_B_LIKE(M,Trajs,[Options])
if (n<4)
  M = args{1};
  trajs = args{2};
  if (n>2)
    Ops = args{3};
  end
  
% LRM_B_LIKE(M,x,Y,Seq,[Options])
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
[Pik, scale] = CalcPik(M,x,Y,Seq,Ops);
Lhood = (sum(log(sum(Pik,2))) + N*log(scale))./prod(size(Y));
[trash, other.C] = max(Pik,[],2);

  
  




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
P = M.order+1;
mlen = max(diff(Seq));
Pik = zeros(n,K);
Piid = zeros(N,D);
Piimk = zeros(N,NumSamps,K);
S = M.S;
TotalSamps = 0;
tries = 1;

while (1)
  TotalSamps = TotalSamps + NumSamps;
  Pik(:) = 0;

  % calculate the density at sampled points
  for k=1:K
    b = randn(NumSamps,1).*sqrt(S(k));
    for j=1:NumSamps
      Xhat = regmat(x-b(j),P-1);
      for d=1:D
        Piid(:,d) = normpdf(Y(:,d),Xhat*M.Mu(:,d,k),M.Sigma(k,d));
      end
      Piimk(:,j,k) = prod(Piid,2);  % get rid of D here (mult it out)
    end
  end
  
  % now scale the data to avoid underflow with long curves
  % and sum across the sample integration points
  scale = mean(mean(mean(Piimk)));
  Piimk_scl = Piimk./scale;  % we don't scale across D; we got rid of D above.
  for k=1:K
    for j=1:TotalSamps
      Pik(:,k) = Pik(:,k) + sprod(Piimk_scl(:,j,k),Seq,mlen);
    end
  end
  clear Piimk_scl;
  Pik = (Pik./TotalSamps) .* (ones(n,1)*M.Alpha');
  if (all(sum(Pik,2))), break; end  % check for stable integration

  % we have detected some zeros, try again?
  if (tries==MaxTries)
    fprintf(['lrm_tt_sh_like: Integration failed, using realmin*1e100 ',...
      'instead.\n']);
    zero = find(sum(Pik,2)==0);
    Pik(zero,:) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
    break;
  else
    fprintf(['lrm_tt_sh_like: Zero membership detected, trying ', ...
        'integration again: %d\n'],tries);
    tries = tries+1;
    S = VarScale*S;  % biased, but gets over some tricky integrations
    NumSamps = floor(Scale*NumSamps);
    Piimk = [zeros(N,NumSamps,K) Piimk]; % save current values
  end
end
