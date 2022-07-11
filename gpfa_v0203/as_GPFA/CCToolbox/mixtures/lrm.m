function M = lrm(trajs,K,order,Ops)
%LRM  Cluster curves w/ multivariate (output) linear regression mixture.
%
%   Model = LRM(Trajs,K,order,[Options])
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find
%    - order : order of polynomial regression
%
%   DefOps = LRM('options) returns the default options.

% Scott Gaffney   29 October 1998
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'lrm';
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
M.zero  = Ops.zero;
M.method = METHOD;

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
N = Seq(end)-1;
% make index for 3-D diagonal variance entries
D = size(Y,2);
diagi = (1 + (D*(0:D-1)) + (0:D-1))' * ones(1,M.K);
diagi = diagi + ones(D,1) * (D^2*(0:M.K-1));

% Initialize the algorithm
M = InitE(M,X,Y,Seq);
M = Mstep(M,X,Y,Seq);


%%%%%%%%%%%%%%%%%%% E-M Algorithm
while (1)
  NumIter = NumIter + 1;
  if(Ops.MsgHnd>=0)
    msgbar(Ops.MsgHnd,sprintf('%sIteration %d',Ops.MsgPrefix,NumIter));
  end
  if (Ops.ShowGraphics), M.Options = showmodel(M,trajs); end

  % Enforce a minimum variance value
  badi = diagi(find(M.Sigma(diagi)<Ops.minvar));
  if (~isempty(badi))
      %fprintf([PROGNAME,': an illegal variance was detected\n']);
      M.Sigma(badi) = Ops.minvar;
  end

  %%% E-Step
  M = Estep(M,X,Y,Seq);
  [Lhood(NumIter),M] = CalcLike(M,N,PROGNAME);
  if (StoppingCondition(Lhood,NumIter,Ops.IterLimit,Ops.stopval))
    break;
  end

  %%% M-Step
  M = Mstep(M,X,Y,Seq);
end
%%%%%%%%%%%%%%%%%%% E-M Algorithm


M = permuteModel(M);
M.Lhood = Lhood(1:NumIter);
[trash, M.C] = max(M.Pik,[],2);
M.NumPoints = prod(size(Y));
M.TrainLhood = M.Lhood(end);
M.TrainLhood_ppt = M.Lhood(end)./M.NumPoints;

% Calculate number of independent parameters
[P,K,D] = size(M.Mu);
if (M.Options.Sigma.Diagonal==1)
  M.NumIndParams = (K-1) + K*P*D + K*D;   % alpha, mu, sigma
else
  M.NumIndParams = (K-1) + K*P*D + K*D*(D+1)/2;
end

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
[N,D] = size(Y);
n = length(Seq)-1;
mlen = max(diff(Seq));

Piik = zeros(N,K);
for k=1:K
  Piik(:,k) = mvnormpdf(Y,X*M.Mu(:,:,k),M.Sigma(:,:,k));
end

M.scale = mean(mean(Piik));
Piik = Piik ./ M.scale;
for k=1:K
  M.Pik(:,k) = sprod(Piik(:,k),Seq,mlen);
end
M.Pik = M.Pik .* (ones(n,1)*M.Alpha');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcLike
%%
function [Lhood,M] = CalcLike(M,N,PROGNAME)
[n,K] = size(M.Pik);
s = sum(M.Pik,2);
if (~all(s))
  fprintf([PROGNAME, ': log(0) detected, using log(K*realmin*1e100).\n']);
  z = find(s==0);
  Pik(z,:) = realmin*1e100*(ones(length(z),1)*M.Alpha');
  s(z) = sum(Pik(z,:),2);
end
Lhood = sum(log(s)) + N*log(M.scale);
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

if (NumIter>1)
  if (isnan(Lhood(NumIter)))
    fprintf('the log-likelihood is equal to NaN.\n');
    dostop = 1;
  end
  if (Lhood(NumIter) < Lhood(NumIter-1))
    %fprintf(['lrm: the log-likelihood appears to have decreased', ...
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
[n,K] = size(M.Pik);
M.Alpha = sum(M.Pik)'./ n;
lens = diff(Seq);

for k=1:K
  [M.Mu(:,:,k), M.Sigma(:,:,k)] = wls(X,Y,copy(M.Pik(:,k),lens));
  if (M.Options.Sigma.Diagonal)
    M.Sigma(:,:,k) = diag(diag(M.Sigma(:,:,k)));
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% InitE
%%
function M = InitE(M,X,Y,Seq)
n = length(Seq)-1;
K = M.K;
Pik = exprnd(.5,n,K);
Pik = Pik ./ (sum(Pik,2)*ones(1,K));
M.Pik = Pik;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%

function Ops = DefaultOptions(Ops);
Ops = SetFieldDef(Ops,'zero','nozero');
Ops = SetFieldDef(Ops,'stopval',1e-5);
Ops = SetFieldDef(Ops,'minvar',1e-5);
Ops = SetFieldDef(Ops,'IterLimit',50);
Ops = SetFieldDef(Ops,'NumEMStarts',1);
Ops = SetFieldDef(Ops,'MsgPrefix','');
Ops = SetFieldDef(Ops,'ShowGraphics',0);
Ops = SetFieldDef(Ops,'MsgHnd',[]);
Ops = SetFieldDef(Ops,'MinLen',[]);
Ops = SetFieldDef(Ops,'Sigma',[]);
Ops.Sigma = SetFieldDef(Ops.Sigma,'Share',0);
Ops.Sigma = SetFieldDef(Ops.Sigma,'Diagonal',1);

function Ops = showmodel(M,trajs)
M = permuteModel(M);
[trash, M.C] = max(M.Pik,[],2);
Ops = viewmodel(M,trajs,M.Options);

function M = permuteModel(M)
M.Mu  = permute(M.Mu,[1 3 2]);
