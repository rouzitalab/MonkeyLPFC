function [Lhood,other] = lrm_like(varargin)
%LRM_LIKE  Calculate log-likelihood with LRM model.
%
%   [Lhood,Other] = LRM_LIKE(M,Trajs,[Options])
%    - M       : trained model
%    - Trajs   : 'Trajs' structure; See also CCTOOLBOX
%    - Options : see MODEL_LIKE
%
%   [Lhood,Other] = LRM_LIKE(M,X,Y,Seq,[Options])
%    - M       : trained model
%    - X,Y,Seq : curves in Sequence format; See also CCTOOLBOX
%    - Options : see MODEL_LIKE
%
%   Other
%     .C   : classification labels
%     .Pik : membership probabilities

% Scott Gaffney   15 February 1999
% DataLab@UCI
% Department of Information and Computer Science
% University of California, Irvine, USA.

PROGNAME = 'lrm_like';
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
% LRM_LIKE(M,Trajs,[Options])
if (n<4)
  M = args{1};
  trajs = args{2};
  if (n>2)
    Ops = args{3};
  end
  
% LRM_LIKE(M,X,Y,Seq,[Options])
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
mlen = max(diff(Seq));
Mu = permute(M.Mu,[1 3 2]);
Piik = zeros(N,K);

% Calculate probability of the data
for j=1:K
  Piik(:,j) = mvnormpdf(Y,X*Mu(:,:,j),M.Sigma(:,:,j));
end

%% Scale the data
scale = mean(mean(Piik));
Piik = Piik./scale;

% calc the loglike
for k=1:K
  Piik(1:n,k) = sprod(Piik(:,k),Seq,mlen);
end
Piik(1:n,:) = Piik(1:n,:) .* (ones(n,1)*M.Alpha');
s = sum(Piik(1:n,:),2);
zero = find(s==0);
if (~isempty(zero))
  Piik(zero,:) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
  s(zero) = sum(Piik(zero,:),2);
end
Lhood = sum(log(s)) + N*log(scale);
Lhood = Lhood./NumPoints;

% Classify the sequences
if (nargout>1)
  [trash, other.C] = max(Piik(1:n,:),[],2);
  other.Pik = Piik(1:n,:)./(s*ones(1,K));
end
