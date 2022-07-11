function M = lrm_d(trajs,K,order,Ops)
%LRM_D  LRM with transformation y=[x]+d
%
%   Model = LRM_D(Trajs,K,order,[Options])
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find
%    - order : order of polynomial regression
%
%   DefOps = LRM_D('options) returns the default options.

% Scott Gaffney   8 May 2003
% Department of Information and Computer Science
% University of California, Irvine
%
% Changes
% -----------------

PROGNAME = 'lrm_d';
METHOD = PROGNAME;
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


%%% Handle Argument Processing
%%
n = nargin;
%
% Check for calling convention
%
% ('options')
%
if (n==1 & strcmp(trajs,'options'))
  M = DefaultOptions([]);
  return;
end

% (Trajs,K,order,[Ops])
if (n<3), error([PROGNAME,': incorrect number of arguments.']); end
M.K = K;
M.order = order;
Ops = cexist('Ops',[]);
%%
%%% End Argument Processing


% Preprocessing
Ops = DefaultOptions(Ops);
M.Options = Ops;
M.method = METHOD;
M.zero = Ops.zero;

% Build data matrices
[Y,X,Seq] = trajs2seq(trajs,M.zero,M.Options.MinLen);
if (size(X,2)~=M.order+1)
  X = regmat(X,M.order);
end

%% Handle graphics output
if (~isfield(Ops,'MsgHnd'))
  Ops.MsgHnd = -1;
end
CreatedMsgBar=0;
if (isempty(Ops.MsgHnd))
  Ops.MsgHnd = msgbar([],'');
  CreatedMsgBar=1;
end
if (Ops.ShowGraphics & isempty(trajs))
  trajs = seq2cell(Y,Seq);
elseif (~Ops.ShowGraphics)
  clear trajs;
end



%***************************************************************************
%   Begin Main Function
%***************************************************************************

%% Define some stuff
NumIter = 0;
Lhood = zeros(1,Ops.IterLimit);

% Initialize the algorithm
M = InitE(M,X,Y,Seq);
M = InitM(M,X,Y,Seq);


%%%%%%%%%%%%%%%%%%% E-M Algorithm
while (1)
  NumIter = NumIter + 1;
  if(Ops.MsgHnd>=0)
    msgbar(Ops.MsgHnd,sprintf('%sIteration %d',Ops.MsgPrefix,NumIter));
  end
  if (Ops.ShowGraphics), M.Options = showmodel(M,trajs); end

  %%% E-Step
  M = Estep(M,X,Y,Seq);
  [Lhood(NumIter), M] = CalcLike(M,PROGNAME);
  if (StoppingCondition(Lhood,NumIter,Ops.IterLimit,Ops.stopval))
    break;
  end

  %%% M-Step
  M = Mstep(M,X,Y,Seq);
end
%%%%%%%%%%%%%%%%%%% E-M Algorithm


M = permuteModel(M); 
M.Lhood = Lhood(1:NumIter);
M.NumPoints = prod(size(Y));
[trash, M.C] = max(M.Pik,[],2);
M.TrainLhood = M.Lhood(end);
M.TrainLhood_ppt = M.Lhood(end)./M.NumPoints;

% Calculate number of independent parameters
[P,K,D] = size(M.Mu);
M.NumIndParams = (K-1) + K*P*D + 2*K*D;  % alpha, mu, sigma, s

if (CreatedMsgBar)
  delete(Ops.MsgHnd);
end


%***************************************************************************
%   End Main Function
%***************************************************************************




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estep
%%
function M = Estep(M,X,Y,Seq)
[P,D,K] = size(M.Mu);
n       = length(Seq)-1;

Pik = zeros(n,K,D);
for k=1:K
  for d=1:D
    Mu       = M.Mu(:,d,k);
    sigma    = M.Sigma(k,d);
    s        = M.S(k,d);
    for i=1:n
      indx   = Seq(i):Seq(i+1)-1;
      ni    = length(indx);
      XMu    = X(indx,:)*Mu;
      M.Vf(i,d,k)  = s*sigma/(ni*s+sigma);
      M.Ef(i,d,k)  = s/(ni*s+sigma)*sum(Y(indx,d)-XMu);
      %M.Ef(i,d,k)  = M.Vf(i,d,k)*sum(Y(indx,d)-XMu)./sigma;
      
      % Pik
      iS = eye(ni)/sigma - 1/(ni*sigma + sigma^2/s);
      Pik(i,k,d) = mvnormpdf_inv(Y(indx,d)',XMu',iS);
    end
  end
end
M.Pik = prod(Pik,3) .* (ones(n,1)*M.Alpha');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcLike
%%
function [Lhood, M] = CalcLike(M,PROGNAME)
[n,K] = size(M.Pik);
s = sum(M.Pik,2);
if (~all(s))
  fprintf([PROGNAME, ': log(0) detected, using log(K*realmin*1e100).\n']);
  zero = find(s==0);
  M.Pik(zero,:) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
  s(zero) = sum(M.Pik(zero,:),2);
end
Lhood = sum(log(s));
M.Pik = M.Pik ./ (s*ones(1,K));  % normalize the memberships while were at it




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% StoppingCondition
%%
function dostop = StoppingCondition(Lhood,NumIter,Limit,Stopval)
dostop=0;

if (NumIter >= Limit)
  dostop = 1;
  return;
end

if (NumIter ~=1)
  if (isnan(Lhood(NumIter)))
    fprintf('the log-likelihood is equal to NaN.\n');
    dostop = 1;
  end
  if (Lhood(NumIter) < Lhood(NumIter-1))
    %fprintf(['the log-likelihood appears to have decreased', ...
    %  ' on this iteration.\n']);
    dostop = 1;
  else
    abs_change = Lhood(NumIter)-Lhood(1);
    if (abs_change==0)
      dostop = 1;
    else
      delta = (Lhood(NumIter)-Lhood(NumIter-1)) / abs_change;
      if (abs(delta) < Stopval)
        dostop = 1;
      end
    end
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mstep
%%
function M = Mstep(M,X,Y,Seq)
[N,D] = size(Y);
K = M.K;
P = M.order+1;
n = length(Seq)-1;
lens = diff(Seq)';

Wk = sum(M.Pik)';

%% Alpha
M.Alpha = Wk./n;

for k=1:K
  Piik = copy(M.Pik(:,k),lens);
  PikX = Piik*ones(1,P).*X;
  f = copy(M.Ef(:,:,k), lens);

  %% s
  M.S(k,:) = M.Pik(:,k)'*(M.Ef(:,:,k).^2 + M.Vf(:,:,k));
  M.S(k,:) = M.S(k,:)./Wk(k);
  
  %% Mu
  M.Mu(:,:,k) = (PikX'*X) \ (PikX'*(Y-f));
    
  %% Sigma
  YxMu = Y-X*M.Mu(:,:,k)-f;
  trVeps = (M.Pik(:,k).*lens)'*M.Vf(:,:,k);
  M.Sigma(k,:) = (Piik'*YxMu.^2 + trVeps) ./ sum(M.Pik(:,k).*lens);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% InitM
%%
function M = InitM(M,X,Y,Seq)
M = Mstep(M,X,Y,Seq);

if (M.Options.InitializeWithExamples)
% this is key for this method
[P,D,K] = size(M.Mu);
rnd = randperm(length(Seq)-1);
i=0;
for k=1:K
  for d=1:D
    i=i+1;
    indx = Seq(rnd(i)):Seq(rnd(i)+1)-1;
    M.Mu(:,d,k) = wls(X(indx,:),Y(indx,d));
  end
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% InitE
%%
function M = InitE(M,X,Y,Seq)
[N,D] = size(Y);
[N,P] = size(X);
n = length(Seq)-1;
K = M.K;

% E-step vars
M.Ef   = zeros(n,D,K);
M.Vf   = zeros(n,D,K);

M.Ef  = randn(n,D,K)*.5;
M.Vf  = rand(n,D,K)*2;
M.Pik = exprnd(.5,n,K);
M.Pik = M.Pik ./ (sum(M.Pik,2)*ones(1,K));






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%
% Options Handling

function Ops = DefaultOptions(Ops);
Ops = SetFieldDef(Ops,'zero','nozero');
Ops = SetFieldDef(Ops,'stopval',1e-5);
Ops = SetFieldDef(Ops,'IterLimit',50);
Ops = SetFieldDef(Ops,'NumEMStarts',1);
Ops = SetFieldDef(Ops,'MsgPrefix','');
Ops = SetFieldDef(Ops,'ShowGraphics',0);
Ops = SetFieldDef(Ops,'MsgHnd',[]);
Ops = SetFieldDef(Ops,'MinLen',[]);
Ops = SetFieldDef(Ops,'InitializeWithExamples',1);

function Ops = showmodel(M,trajs)
M.Mu = permute(M.Mu, [1 3 2]);  % make p-K-D
M.Ef  = permute(M.Ef,[1 3 2]);
[trash, M.C] = max(M.Pik,[],2);
Ops = viewmodel(M,trajs,M.Options);

function M = permuteModel(M)
M.Ef  = permute(M.Ef,[1 3 2]);
M.Vf  = permute(M.Vf,[1 3 2]);
M.Mu  = permute(M.Mu,[1 3 2]);

