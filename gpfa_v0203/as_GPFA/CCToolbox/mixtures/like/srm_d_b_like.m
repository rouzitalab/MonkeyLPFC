function [Lhood,other] = srm_d_b_like(varargin)
%SRM_D_B_LIKE  Calculate log-likelihood with SRM_D_B model.
%
%   [Lhood,Other] = SRM_D_B_LIKE(M,Trajs,[Options])
%    - M       : trained model
%    - Trajs   : 'Trajs' structure; See also CCTOOLBOX
%    - Options : see MODEL_LIKE
%
%   [Lhood,Other] = SRM_D_B_LIKE(M,X,Y,Seq,[Options])
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

PROGNAME = 'srm_d_b_like';
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
% SRM_D_B_LIKE(M,Trajs,[Options])
if (n<4)
  M = args{1};
  trajs = args{2};
  if (n>2)
    Ops = args{3};
  end
  
% SRM_D_B_LIKE(M,x,Y,Seq,[Options])
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
  Ops.NumSamps = 300;
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
function Pik = CalcPik(M,x,Y,Seq,Ops)
% Numerical integration
NumSamps = Ops.NumSamps;
MaxTries = Ops.MaxAtts;
Scale = Ops.Scale;
VarScale = Ops.VarScale;

[N,D] = size(Y);
n = length(Seq)-1;
K = M.K;
Pid = zeros(n,D);
Pik = zeros(n,K);
S = M.S;  T = M.T;
TotalSamps = 0;
tries = 1;

while (1)
  TotalSamps = TotalSamps + NumSamps;

  % calculate the density at sampled points
  for k=1:K
    b = randn(NumSamps,1).*sqrt(S(k));      % sample from N(0,s)
    for j=1:NumSamps
      Xhat = bsplinebasis(M.knots,M.order,x-b(j));
      for d=1:D
        sigma = M.Sigma(k,d);  t = T(k,d);
        for i=1:n
          indx = Seq(i):Seq(i+1)-1;  ni = length(indx);
          iV = eye(ni)./sigma - 1/(ni*sigma + sigma^2/t);
          Pid(i,d) = mvnormpdf_inv(Y(indx,d)',Xhat(indx,:)*M.Mu(:,d,k),iV);
        end
      end
      Pik(:,k) = Pik(:,k) + prod(Pid,2);
    end
  end
  pik = (Pik./TotalSamps) .* (ones(n,1)*M.Alpha');  % we keep Pik for next try!
  if (all(sum(pik,2))), break; end

  % we have detected some zeros, try again?
  if (tries==MaxTries)
    fprintf('Integration failed, using realmin*1e100 instead.\n');
    zero = find(sum(pik,2)==0);
    pik(zero,:) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
    break;
  else
    fprintf('Zero membership detected, trying integration again: %d\n',tries);
    tries = tries+1;
    NumSamps = round(NumSamps*Scale);
    S = VarScale*S;  % biased, but gets over some tricky integrations
    %T = 1.25*T;
  end
end
