function ModelOut = rerm(X,Y,K,Seq,Options)
%RERM  Learn a Random Effects Regression Mixture model
%
%   IMPORTANT: This is old experimental code. The algorithm works fine
%   but there may be problems with miscellaneous pieces of code, e.g.,
%   graphics output, error checking, etc.
%
%   ModelOut = RERM(X,Y,K,SEQ,[OPTIONS]) 
%   clusters the curve data in X and Y by forming a mixture density with K 
%   components, where each component is a linear regression model.
%
%   SEQ is an index vector which gives the positions in X 
%   where different sequences of points exist. For example, 
%   SEQ = [1 10 23 56] states that the first sequence consists
%   of 9 points located at X(1:10-1), the second sequence 
%   consists of 13 points located at X(10:23-1), and the third
%   and last sequence consists of 33 points located at
%   X(23:56-1).  Note that LENGTH(SEQ) = (1 + number of sequences).
%   If your data does not contain sequences, or if you do not
%   wish to use this information (or do not know this information
%   apriori), then simply set SEQ = [].
%
%   RERM outputs a model structure. See rerm_model() for more info.
%

% Scott Gaffney   25 July 2002
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'rerm';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

if (nargin < 3)
  error([PROGNAME ': incorrect number of arguments passed to this function.']);
end

[xN,P] = size(X);
[yN,D] = size(Y);
if (xN~=yN)
  error([PROGNAME ': the number of rows in X must equal those in Y.']);
end
if (exist('Seq','var') ~= 1)
  Seq = [1:yN yN+1];  % no sequence information
end

% Options Handling
Options = cexist('Options',[]);
Options = HandleOptions(Options,P);





%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*
%* Begin Main Function
%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*

%% Define some stuff
NumIter = 0;
ITER_LIMIT = 100;
Lhood = zeros(1,ITER_LIMIT);

% Initialize the algorithm
[Sigma,Alpha,Mu,R] = InitEM(Options,X,Y,K,Seq);


%%%%%%%%%%%%%%%%%%% E-M Algorithm
while (1)
  NumIter = NumIter + 1;

  %%% E-Step
  [Pik,Beta_ikd,Vikd] = Estep(Options,X,Y,Seq,Sigma,Alpha,Mu,R);
  [Lhood(NumIter),Pik] = CalcLike(Options,Pik);
  if (StoppingCondition(Lhood,NumIter,ITER_LIMIT))
    break;
  end

  %%% M-Step
  [Alpha,Mu,R] = TopLevel_Mstep(Options,Pik,Beta_ikd,Vikd,R);
  [Sigma] = BottomLevel_Mstep(Options,X,Y,Seq,Pik,Beta_ikd,Vikd);
  
  if (0)
    Beta_i = PredictBeta(Pik,Beta_ikd);
    Options=HandleGraphics(Options,X,Y,Seq,Beta_i,Sigma,Pik,Mu,R,NumIter,Lhood);
  elseif(Options.MsgHnd~=-1)
    msgbar(Options.MsgHnd,sprintf('%sIteration %d',Options.MsgPrefix,NumIter));    
  end
end

%%%%%%%%%%%%%%%%%%% E-M Algorithm

Lhood = Lhood(1:NumIter);
[trash, C] = max(Pik,[],2);
NumPoints = prod(size(Y));
Beta_i = PredictBeta(Pik,Beta_ikd);
ModelOut = rerm_model(Beta_i,Sigma,Alpha,Mu,R,Pik,Lhood,C,NumPoints, ...
    Options.back,Options.order,Options.zero,Options.method);
ModelOut.Beta_ikd = Beta_ikd;
ModelOut.TrainLhood = Lhood(end);
ModelOut.TrainLhood_ppt = Lhood(end)/ModelOut.NumPoints;

% GUI Cleanup
if (Options.DoGraphics)
  ShutdownGraphics(Options);
end


%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*
%*  End Main Function
%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*



%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*
%* EM Functions
%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estep (Non-normalized)
%%
function [Pik,Beta_ikd,Vikd] = Estep(Options,X,Y,Seq,Sigma,Alpha,Mu,R)
n = length(Seq)-1;
[N,D] = size(Y);
[P,K,D] = size(Mu);

Pik = zeros(n,K,D);
Beta_ikd = zeros(P,n,K,D);
Vikd = zeros(P,P,n,K,D);
for k=1:K
  for d=1:D
    Rinv = inv(R(:,:,k,d));
    for i=1:n
      indx = Seq(i):Seq(i+1)-1;
      n_i = length(indx);
      V_y = X(indx,:)*R(:,:,k,d)*X(indx,:)' + Sigma(i,d)*eye(n_i);
      Pik(i,k,d) = mvnormpdf(Y(indx,d)',X(indx,:)*Mu(:,k,d),V_y);
      
      A = 1/Sigma(i,d)*X(indx,:)'*X(indx,:) + Rinv;
      c = 1/Sigma(i,d)*X(indx,:)'*Y(indx,d) + Rinv*Mu(:,k,d);
      Beta_ikd(:,i,k,d) = A\c;
      
      Vikd(:,:,i,k,d) = inv(1/Sigma(i,d)*X(indx,:)'*X(indx,:) + Rinv);
    end
  end
end
Pik = prod(Pik,3);  % scaling problems?
Pik = Pik .* (ones(n,1)*Alpha');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BottomLevel_Mstep
%%
function [Sigma] = BottomLevel_Mstep(Options,X,Y,Seq,Pik,Beta_ikd,Vikd)
[n,K] = size(Pik);
[N,D] = size(Y);

Sigma = zeros(n,D);
for i=1:n
  indx = Seq(i):Seq(i+1)-1;
  n_i = length(indx);
  for k=1:K
    for d=1:D
      yxb = Y(indx,d)- X(indx,:)*Beta_ikd(:,i,k,d);
      Sigma(i,d) = Sigma(i,d) + ...
        Pik(i,k)* (yxb'*yxb + trace(X(indx,:)'*X(indx,:)*Vikd(:,:,i,k,d)));
    end
  end
  if (~Options.Sigma.Single)
    Sigma(i,:) = Sigma(i,:)./n_i;
  end
end
if (Options.Sigma.Single)
  Sigma = ones(n,1)*(sum(Sigma)/N);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TopLevel_Mstep
%%
function [Alpha,Mu,R] = TopLevel_Mstep(Options,Pik,Beta_ikd,Vikd,R)
[P,n,K,D] = size(Beta_ikd);

if (Options.back)
  bK=K-1;   % do not learn background parameters Mu, R
else
  bK=K;
end

% Mixing weights (w/ Dirichlet prior)
Alpha = (sum(Pik)'+ Options.Alpha.Dirich_a-1)./ (n+sum(Options.Alpha.Dirich_a)-K);

% Means; do not learn background mean.
for k=1:bK
  for d=1:D
    Mu(:,k,d) = Beta_ikd(:,:,k,d)*Pik(:,k);
    Mu(:,k,d) = Mu(:,k,d) ./ ((ones(P,1)*sum(Pik(:,k)))+ realmin);
  end
end

% Covariance Matrices
Rtmp = zeros(P,P,K,D);
for k=1:K
  P_i = reshape(repmat(Pik(:,k),P,P),n,P,P);
  P_i = permute(P_i,[2 3 1]);
  for d=1:D
    BMu = Beta_ikd(:,:,k,d)-(Mu(:,k,d)*ones(1,n));
    Rtmp(:,:,k,d) = BMu *spdiags(Pik(:,k),0,n,n)* BMu' + sum(Vikd(:,:,:,k,d).*P_i,3);
    
    if (~Options.R.Share & k<=bK)  % do not learn the background
      R(:,:,k,d) = (Rtmp(:,:,k,d) + Options.R.Wish_iB) ./ ...
        (sum(Pik(:,k))+ (Options.R.Wish_a-(P+1)));
      if (Options.R.Diagonal)
        R(:,:,k,d) = diag(diag(R(:,:,k,d)));
      end
    end
  end
end

if (Options.R.Share)
  Rshare = squeeze(sum(Rtmp,3));  clear Rtmp;  % note: we use the background here
  for d=1:D
    Rshare(:,:,d) = (Rshare(:,:,d) + Options.R.Wish_iB) ...
      ./ (n + Options.R.Wish_a-(P+1));
    if (Options.R.Diagonal)
      Rshare(:,:,d) = diag(diag(Rshare(:,:,d)));
    end
    R(:,:,1:bK,d) = repmat(Rshare(:,:,d),[1 1 bK]);  % do not learn the background
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcLike
%%
function [Lhood,Pik] = CalcLike(Options,Pik)
s = sum(Pik,2);
s(find(s==0)) = realmin;
Lhood = sum(log(s));
Pik = Pik ./ (s*ones(1,size(Pik,2)));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% StoppingCondition
%%
function dostop = StoppingCondition(Lhood,NumIter,Limit)
stopval = 1e-7;
dostop=0;

if (NumIter > Limit)
  dostop = 1;
  return;
end

if (NumIter > 3)
  if (isnan(Lhood(NumIter)))
    fprintf('the log-likelihood is equal to NaN.\n\n');
    Lhood(isnan(Lhood)) = -inf;
    dostop = 1;
  elseif (Lhood(NumIter) < Lhood(NumIter-1))
    fprintf('log-likelihood decreased from %g to %g on iteration %d.\n', ...
        Lhood(NumIter-1),Lhood(NumIter),NumIter);
    dostop = 1;
  else
    delta = (Lhood(NumIter)-Lhood(NumIter-1)) / (Lhood(NumIter)-Lhood(2));
    if (abs(delta) < stopval)
      dostop = 1;
    end
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PredictBeta
%%
function Beta_i = PredictBeta(Pik,Beta_ikd)
[P,n,K,D] = size(Beta_ikd);
Pik_PnKD = permute(repmat(Pik,[1 1 P D]),[3 1 2 4]);
Beta_i = squeeze(sum(Beta_ikd.*Pik_PnKD,3));










%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*
%* Initialization Functions
%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% InitEM
%%
function [Sigma,Alpha,Mu,R] = InitEM(Options,X,Y,K,Seq)
if (Options.back | Options.BottomInit)
  [B,Sigma] = InitBottomLevel(Options,X,Y,Seq);
  [Alpha,Mu,R] = InitTopLevel(Options,B,K);
else
  [Sigma,Alpha,Mu,R] = WLSInit(Options,X,Y,K,Seq);
  if (Options.RandomInit)
    Mu = [36 -.9 0; 0 .9 0; 40 -2 .05]';
  end
end

  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WLSInit
%%
function [Sigma,Alpha,Mu,R] = WLSInit(Options,X,Y,K,Seq)
n = length(Seq)-1;
D = size(Y,2);
P = size(X,2);
Pik = RandPik(Options,n,K);
Alpha = ones(K,1)./K;
Mu = zeros(P,K,D);
for k=1:K
  [Mu(:,k,:), XYVar(:,:,k)] = wls(X,Y,copy(Pik(:,k),diff(Seq)));
end
DiagVar = diag3(XYVar);
R = RandR(DiagVar,D,K,P);
Sigma = RandSigma(DiagVar,D,n);
if (Options.Sigma.Single)
  Sigma = (ones(size(Sigma,1),1)*sum(Sigma))/size(Y,1);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RandPik
%%
function [Pik] = RandPik(Options,N,K)
if (~isempty(Options) & Options.back)
  Pik = exprnd(.5,N,K-1);
  Pik = Pik ./ (sum(Pik,2)*ones(1,K-1));
  Pik(:,K) = ones(N,1);
else
  Pik = exprnd(.5,N,K);
  Pik = Pik ./ (sum(Pik,2)*ones(1,K));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RandR
%%
function [R] = RandR(DiagVar,D,K,P)
R = zeros(P,P,K,D);
for d=1:D
  R(:,:,1,d) = diag(diag(randn(P,P) + mean(DiagVar(d,:))));
  for k=2:K
    R(:,:,k,d) = R(:,:,1,d);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RandSigma
%%
function [Sigma] = RandSigma(DiagVar,D,N)
Sigma = randn(N,D) + ones(N,1)*mean(DiagVar');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HandleOptions
%%
function Options = HandleOptions(Options,P)
Options = SetFieldDef(Options,'back',0);
Options = SetFieldDef(Options,'order',[]);
Options = SetFieldDef(Options,'zero',[]);
Options = SetFieldDef(Options,'method','rerm');
Options = SetFieldDef(Options,'R',mkR(P));
Options = SetFieldDef(Options,'WLSInit',1);
Options = SetFieldDef(Options,'BottomInit',0);
Options = SetFieldDef(Options,'MsgHnd',0);
Options = SetFieldDef(Options,'Sigma',[]);
Options.Sigma = SetFieldDef(Options.Sigma,'Single',1);







