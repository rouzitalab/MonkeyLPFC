function M = srm_d(trajs,K,knots,order,Ops)
%SRM_D  SRM with transformation y=[x]+d
%
%   Model = SRM_D(Trajs,K,[knots],order,[Options])
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find
%    - knots : vector of knot values (pass in [] if you want the algorithm
%              to select the knots "automatically"; Use Options.Kn if
%              you would like to specify the number of knots); note that
%              the resulting augmented knot sequence may have a longer
%              length (this is normal).
%    - order : order of polynomial regression
%
%   DefOps = SRM_D('options) returns the default options.
%
%   Options (some fields)
%    .KnotMult : (default .25) multiplied by the complete length of data 
%                interval to get the number of knots (minimum of 1)

% Scott Gaffney   8 May 2003
% Department of Information and Computer Science
% University of California, Irvine
%
% Changes
% -----------------

PROGNAME = 'srm_d';
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
M.method = METHOD;
M.zero = Ops.zero;

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
  [DoStop,Ops] = StoppingCondition(Lhood,NumIter,Ops);
  if (DoStop), break; end

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
      M.Vf(i,d,k)  = sigma*s/(ni*s + sigma);
      M.Ef(i,d,k)  = s/(ni*s + sigma)*sum(Y(indx,d)-XMu);
      
      % Pik
      iS = eye(ni)./sigma - 1/(ni*sigma + sigma^2/s);
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
M.Pik = M.Pik ./ (s*ones(1,K));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% StoppingCondition
%%
function [dostop,Ops] = StoppingCondition(Lhood,NumIter,Ops)
dostop=0;

if (NumIter >= Ops.IterLimit)
  dostop = 1;
  return;
end

if (NumIter ~=1)
  if (isnan(Lhood(NumIter)))
    fprintf('the log-likelihood is equal to NaN.\n');
    dostop = 1;
  end
  if (Lhood(NumIter) < Lhood(NumIter-1))
    Ops.NumDec = Ops.NumDec+1;
    %fprintf(['the log-likelihood wabbled down', ...
    % ' on this iteration from ',num2str(Lhood(NumIter-1)),' to ', ...
    % num2str(Lhood(NumIter)),'.\n']);
    if (Ops.NumDec>=Ops.MaxDec)
      dostop = 1;
    end
  else
    abs_change = Lhood(NumIter)-Lhood(1);
    if (abs_change==0)
      dostop = 1;
    else
      delta = (Lhood(NumIter)-Lhood(NumIter-1)) / abs_change;
      if (abs(delta) < Ops.stopval)
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
P = length(M.knots)-M.order;
n = length(Seq)-1;
lens = diff(Seq)';
mlen = max(lens);

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

  % check valid knot interval
  mass = sum(PikX);
  invalid = find(mass<1);  % simple way to ignore coefs with no data
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
  PikX(:,invalid) = [];
  
  %% Mu
  M.Mu(:,:,k) = (PikX'*X) \ (PikX'*(Y-f));
  
  % fill-in invalid coefficients w/ closest valid ones
  fd = find(invalid<valid(1));   
  M.Mu(invalid(fd),:,k) = ones(length(fd),1)*M.Mu(valid(1),:,k);
  fd = find(invalid>valid(end)); 
  M.Mu(invalid(fd),:,k) = ones(length(fd),1)*M.Mu(valid(end),:,k);
  
  %% Sigma
  YxMu = Y-X*M.Mu(:,:,k)-f;
  trVeps = (M.Pik(:,k).*lens)'*M.Vf(:,:,k);
  M.Sigma(k,:) = (Piik'*YxMu.^2 + trVeps) ./ sum(M.Pik(:,k).*lens);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% InitM
%%
function M = InitM(M,X,Y,Seq)
M = Mstep(M,X,Y,Seq);  % default init

if (M.Options.InitializeWithExamples)
% this is key for this method
[P,D,K] = size(M.Mu);
rnd = randperm(length(Seq)-1);
i=0;
for k=1:K
  i=i+1;
  indx = Seq(rnd(i)):Seq(rnd(i)+1)-1;
  x = X(indx,:);
  
  % check valid knot interval
  mass = sum(x);
  invalid = find(mass<1);
  x(:,invalid) = [];
  valid = (1:P)';  valid(invalid)=[];
  for d=1:D
    M.Mu(valid,d,k) = wls(x,Y(indx,d));
  end
  
  % fill-in invalid coefficients w/ closest valid ones
  fd = find(invalid<valid(1));   
  M.Mu(invalid(fd),:,k) = ones(length(fd),1)*M.Mu(valid(1),:,k);
  fd = find(invalid>valid(end)); 
  M.Mu(invalid(fd),:,k) = ones(length(fd),1)*M.Mu(valid(end),:,k);
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
Ops = SetFieldDef(Ops,'IterLimit',75);
Ops = SetFieldDef(Ops,'MaxDec',2);
Ops = SetFieldDef(Ops,'NumDec',0);
Ops = SetFieldDef(Ops,'NumEMStarts',1);
Ops = SetFieldDef(Ops,'MsgPrefix','');
Ops = SetFieldDef(Ops,'ShowGraphics',0);
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

