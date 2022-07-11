function M = srm(trajs,K,knots,order,Ops)
%SRM  Cluster curves w/ Spline Regression Mixtures (SRM)
%
%   Model = SRM(Trajs,K,knots,order,[Options])
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find
%    - knots : vector of knot values (pass in [] if you want the algorithm
%              to select the knots "automatically"; Use Options.Kn if
%              you would like to specify the number of knots); note that
%              the resulting augmented knot sequence may have a longer
%              length (this is normal).
%    - order : order of polynomial regression
%
%   DefOps = SRM('options) returns the default options.
%
%   Options (some fields)
%    .KnotMult   : (default .25) multiplied by the complete length of data 
%                  interval to get the number of knots (minimum of 1)

% Scott Gaffney   9 October 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'srm';
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
%
% (Trajs,K,knots,order,[Options])
if (n<4), error([PROGNAME,': incorrect number of arguments.']); end
M.K = K;
M.knots = knots;
M.order = order;
Ops = cexist('Ops',[]);
%%
%%% End Argument Processing


% Preprocessing
Ops = DefaultOptions(Ops);
M.Options = Ops;
M.zero  = Ops.zero;
M.method = METHOD;

% build data matrices
[Y,X,Seq] = trajs2seq(trajs,M.zero,M.Options.MinLen);
if (isempty(M.knots))
  M.knots = SelectKnots(M,X,Y,Seq);
end
if (size(X,2)==1)  % problem if X is a row-vector
  X = bsplinebasis(M.knots,M.order,X);
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
M.NumPoints = prod(size(Y));
[trash, M.C] = max(M.Pik,[],2);
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
  zero = find(s==0);
  M.Pik(zero,:) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
  s(zero) = sum(Pik(zero,:),2);
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
[n,K] = size(M.Pik);
P = length(M.knots)-M.order;
lens = diff(Seq)';

%% Alpha
Wk = sum(M.Pik);
M.Alpha = Wk'./ n;

for k=1:K
  sqPiik = copy(sqrt(M.Pik(:,k)),lens);
  sqPikX = sqPiik*ones(1,P).*X;
  
  % check valid knot interval
  mass = sum(sqPikX);
  invalid = find(mass<1);
  valid = (1:P)';  valid(invalid)=[];

  % specialized error code: shouldn't be needed under normal circumstances
  if (isempty(valid))  % this only happens when there isn't enough data...
    fprintf('This cluster has no support.\n');  % ...and in this case, you...
    if (P<=2)          % ...shouldn't be fitting this model, but this...
      valid = [1:P]';  % ...section will get rid of crashing issues in any case
      invalid = [];
    else
      valid = [2:P-1]';  % often this will get over the problem, but the...
      invalid = [1 P]';  % ...resulting Mu won't be accurate
    end
  end
  
  % remove unwanted spline coefs (those with no support)
  sqPikX(:,invalid) = [];

  %% Mu
  M.Mu(valid,:,k) = sqPikX \ (sqPiik*ones(1,D).*Y);

  % fill-in invalid coefficients w/ closest valid ones
  f = find(invalid<valid(1));
  M.Mu(invalid(f),:,k) = ones(length(f),1) * M.Mu(valid(1),:,k);
  f = find(invalid>valid(end));
  M.Mu(invalid(f),:,k) = ones(length(f),1) * M.Mu(valid(end),:,k);

  %% Sigma
  Piik = copy(M.Pik(:,k),lens);
  YxMu = Y-X*M.Mu(:,:,k);
  M.Sigma(:,:,k) = ( (Piik*ones(1,D).*YxMu)'*YxMu )./sum(lens.*M.Pik(:,k));

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
M.Pik = exprnd(.5,n,K);
M.Pik = M.Pik ./ (sum(M.Pik,2)*ones(1,K));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SelectKnots
%%
function knots = SelectKnots(M,x,Y,Seq)
Start = min(x(Seq(1:end-1)));
End = max(x(Seq(2:end)-1));
len = max(diff(Seq));

knotmult = .25;
if (isfield(M.Options,'KnotMult'))
  knotmult = M.Options.KnotMult;
end

% check for number of knots
if (isfield(M.Options,'Kn') & ~isempty(M.Options.Kn))
  Kn = M.Options.Kn;
else
  Kn = max(ceil(len*knotmult),1);
end
knots = linspace(Start,End,Kn);  % uniform spacing
knots = addendpts(knots,M.order);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%

function Ops = DefaultOptions(Ops);
Ops = SetFieldDef(Ops,'zero','nozero');
Ops = SetFieldDef(Ops,'stopval',1e-7);
Ops = SetFieldDef(Ops,'IterLimit',50);
Ops = SetFieldDef(Ops,'NumEMStarts',1);
Ops = SetFieldDef(Ops,'MsgPrefix','');
Ops = SetFieldDef(Ops,'ShowGraphics',0);
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
