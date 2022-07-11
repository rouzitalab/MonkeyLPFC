function M = srm_ab(trajs,K,knots,order,Ops)
%SRM_AB  SRM with transformation y=[ax+b]
%
%   Model = SRM_AB(Trajs,K,knots,order,[Options])
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find
%    - knots : vector of knot values (pass in [] if you want the algorithm
%              to select the knots "automatically"; Use Options.Kn if
%              you would like to specify the number of knots); note that
%              the resulting augmented knot sequence may have a longer
%              length (this is normal).
%    - order : order of polynomial regression
%
%   DefOps = SRM_AB('options) returns the default options.
%
%   Options (some fields)
%    .Interval  : range of allowed translations (e.g., [-2 2]), you  must
%                  provide this if knots is empty when this functions is called.
%    .KnotMult : (default .25) multiplied by the complete length of data 
%                interval to get the number of knots (minimum of 1)

% Scott Gaffney   1 Aug 2003
% Department of Information and Computer Science
% University of California, Irvine
%
% Changes
% -----------------

PROGNAME = 'srm_ab';
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
[Y,x,Seq] = trajs2seq(trajs,M.zero,M.Options.MinLen);
if (isempty(M.knots))
  if (isempty(Ops.Interval))
    errorbox(['You must provide either ', ...
      'knots or Options.Interval.'],[PROGNAME,':   Argument Error']);
    return;
  end
  M.knots = SelectKnots(M,x,Y,Seq);
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
M = InitE(M,x,Y,Seq);
M = InitM(M,x,Y,Seq);


%%%%%%%%%%%%%%%%%%% E-M Algorithm
while (1)
  NumIter = NumIter + 1;
  if(Ops.MsgHnd>=0)
    msgbar(Ops.MsgHnd,sprintf('%sIteration %d',Ops.MsgPrefix,NumIter));
  end
  if (Ops.ShowGraphics), M.Options = showmodel(M,trajs); end

  %%% E-Step
  M.NumIter = NumIter;
  M = Estep(M,x,Y,Seq);
  [Lhood(NumIter), M] = CalcLike(M,N);
  [DoStop,Ops] = StoppingCondition(Lhood,NumIter,Ops);
  if (DoStop), break; end

  %%% M-Step
  M = Mstep(M,x,Y,Seq);
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
M.NumIndParams = (K-1) + K*P*D + K*D + 2*K;  % alpha, mu, sigma, s, r
if (M.Options.Sigma.Share==1)
  M.NumIndParams = M.NumIndParams - K*(D-1);
end

if (CreatedMsgBar)
  delete(M.Options.MsgHnd);
end


%***************************************************************************
%   End Main Function
%***************************************************************************




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% postval 
%%
%%
function val = postval(pt,x,Y,Mu,knots,order,r,s,sigma)
a = pt(1);  b = pt(2);
Xhat = bsplinebasis(knots,order,a*x-b);
val = sum(sum((Y-Xhat*Mu).^2)./sigma) + (a-1)^2/r + b^2/s;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estep
%%
% function M = Estep(M,x,Y,Seq,trajs,NumIter)
function M = Estep(M,x,Y,Seq)
[P,D,K] = size(M.Mu);
n = length(Seq)-1;
fun = @postval;
SearchOps = M.Options.SearchOps;

for k=1:K
  r        = M.R(k);
  s        = M.S(k);
  Mu       = M.Mu(:,:,k);
  sigma    = M.Sigma(k,:);
  for i=1:n
    indx   = Seq(i):Seq(i+1)-1;

    % find the mode of the posterior and set Ea and Eb to this point
    a0 = M.Ea(i,k) - (1-M.Options.PropStart)*(M.Ea(i,k)-1);
    b0 = M.Eb(i,k) * M.Options.PropStart;
    pt0 = [a0 b0]';
    maxpt = fminsearch(fun,pt0,SearchOps,x(indx),Y(indx,:),Mu, ...
      M.knots,M.order,r,s,sigma);
    a = maxpt(1);  b = maxpt(2);
    if (a==0), fprintf('a was zero.\n'); a=M.Ea(i,k); end
    M.Ea(i,k) = a;  M.Eb(i,k) = b;

    % Now we will calc the inverse Information at Ea and Eb.
    % First we make Xhat at (Ea,Eb) and then set up the gradients...
    Xhat = bsplinebasis(M.knots,M.order,a*x(indx)-b);
    Dx  = bsplinebasis(M.knots,M.order,a*x(indx)-b,1);
    D2x = bsplinebasis(M.knots,M.order,a*x(indx)-b,2);
    
    % ...then make repeat calculations
    DxMu = Dx*Mu;  % ni x D
    DxMux = DxMu.*x(indx,ones(1,D));  % ni x D
    D2xMu = D2x*Mu;  % ni x D
    D2xMux = D2xMu.*x(indx,ones(1,D));  % ni x D
    YxMu = Y(indx,:)-Xhat*Mu;  % ni x D
    
    % ...and now we can calculate the information.
    Ia  = sum(sum(YxMu.*(D2xMux.*x(indx,ones(1,D))))./sigma) ...
         -sum(sum(DxMux.^2)./sigma) - 1./r;
    Ib  = sum(sum(YxMu.*D2xMu)./sigma) - sum(sum(DxMu.^2)./sigma) - 1./s;
    Iab = -sum(sum(YxMu.*D2xMux)./sigma) + sum(sum(DxMu.*DxMux)./sigma);
    
    % the variance is the inverse of the negative information
    V = inv(-[Ia Iab; Iab Ib]);
    M.Va(i,k)  = V(1);
    M.Vb(i,k)  = V(4);
    M.Vab(i,k) = V(2);
  end
end

if (any(M.Va<=0) | any(M.Vb<=0))
  fprintf('A negative variance was detected....\n');
  fa = find(M.Va<=0);    fb = find(M.Vb<=0);
  M.Va(fa) = -M.Va(fa);  M.Vb(fb) = -M.Vb(fb);
end

M = CalcPik(M,x,Y,Seq);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcPik 
%%
function M = CalcPik(M,x,Y,Seq)
% Numerical integration
NumSamps = 100;
MaxTries = 5;
[N,D] = size(Y);
n = length(Seq)-1;
K = M.K;
mlen = max(diff(Seq));
Piid = zeros(N,D);
Piimk = zeros(N,NumSamps,K);
M.Pik(:) = 0;
S = M.S;  R = M.R;
TotalSamps = 0;
tries = 1;

while (1)
  TotalSamps = TotalSamps + NumSamps;
  M.Pik(:) = 0;

  % calculate the density at sampled points
  for k=1:K
    r = R(k);   s = S(k);
    a = randn(NumSamps,1).*sqrt(r) + 1;  % sample from N(1,r)
    b = randn(NumSamps,1).*sqrt(s);      % sample from N(0,s)
    for j=1:NumSamps
      Xhat = bsplinebasis(M.knots,M.order,a(j)*x-b(j));
      for d=1:D
        Piid(:,d) = normpdf(Y(:,d),Xhat*M.Mu(:,d,k),M.Sigma(k,d));
      end
      Piimk(:,j,k) = prod(Piid,2);
    end
  end
  
  % now scale the data to avoid underflow with long curves
  % and sum across the sample integration points
  M.scale = mean(mean(mean(Piimk)));
  Piimk_scl = Piimk./M.scale;
  for k=1:K
    for j=1:TotalSamps
      M.Pik(:,k) = M.Pik(:,k) + sprod(Piimk_scl(:,j,k),Seq,mlen);
    end
  end
  clear Piimk_scl;
  M.Pik = (M.Pik./TotalSamps) .* (ones(n,1)*M.Alpha');
  if (all(sum(M.Pik,2))), break; end

  % we have detected some zeros, try again?
  if (tries==MaxTries)
    fprintf(['srm_ta_sh: Integration failed, using realmin*1e100 ',...
      'instead.\n']);
    zero = find(sum(M.Pik,2)==0);
    M.Pik(zero,:) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
    break;
  else
    fprintf(['srm_ta_sh: Zero membership detected, trying ', ...
        'integration again: %d\n'],tries);
    tries = tries+1;
    S = 1.25*S;  % biased, but gets over some tricky integrations
    R = 1.25*R;  % biased, but gets over some tricky integrations
    Piimk = [zeros(N,NumSamps,K) Piimk]; % save current values
  end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcLike
%%
function [Lhood, M] = CalcLike(M,N)
[n,K] = size(M.Pik);
s = sum(M.Pik,2);
Lhood = sum(log(s)) + N*log(M.scale);
M.Pik = M.Pik ./ (s*ones(1,K));  % normalize the memberships while were at it




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
    fprintf('the log-likelihood is equal to NaN.');
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
function M = Mstep(M,x,Y,Seq)
P = length(M.knots)-M.order;
K = M.K;
[N,D] = size(Y);
n = length(Seq)-1;
lens = diff(Seq)';
Wk = sum(M.Pik)';

%% Alpha
M.Alpha = Wk./n;

% start the M-loop
for k=1:K
  Piik = copy(M.Pik(:,k),lens);

  %% s
  M.S(k) = M.Pik(:,k)'*(M.Eb(:,k).^2 + M.Vb(:,k));
  M.S(k) = M.S(k)./Wk(k);

  %% r
  M.R(k) = M.Pik(:,k)'*((M.Ea(:,k)-1).^2 + M.Va(:,k));
  M.R(k) = M.R(k)./Wk(k);

  % now make Xhat and finish up the M-Step
  Xhat = bsplinebasis(M.knots,M.order, ...
    copy(M.Ea(:,k),lens).*x- copy(M.Eb(:,k),lens));
  PikXhat = Piik*ones(1,P).*Xhat;   

  % check valid knot interval
  mass = sum(PikXhat);
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
  PikXhat(:,invalid) = [];
  
  %% Mu
  M.Mu(:,:,k) = (PikXhat'*Xhat) \ (PikXhat'*Y);
  
  % fill-in invalid coefficients w/ closest valid ones
  fd = find(invalid<valid(1));
  M.Mu(invalid(fd),:,k) = ones(length(fd),1) * M.Mu(valid(1),:,k);
  fd = find(invalid>valid(end));
  M.Mu(invalid(fd),:,k) = ones(length(fd),1) * M.Mu(valid(end),:,k);

  % Sigma
  M.Sigma(k,:) = Piik'*((Y-Xhat*M.Mu(:,:,k)).^2)./ sum(M.Pik(:,k).*lens);
end
if (M.Options.Sigma.Share)
  M.Sigma = (sum(M.Sigma,2)./D)*ones(1,D);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% InitM
%%
function M = InitM(M,x,Y,Seq)
M = Mstep(M,x,Y,Seq);  % default init

if (M.Options.InitializeWithExamples)
[P,D,K] = size(M.Mu);
X = bsplinebasis(M.knots,M.order,x);
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
function M = InitE(M,x,Y,Seq)
[N,D] = size(Y);
n = length(Seq)-1;
K = M.K;
P = length(M.knots)-M.order;

% E-step vars
% M.Ea   = zeros(n,K);
% M.Eb   = zeros(n,K);
M.Va   = zeros(n,K);
M.Vb   = zeros(n,K);
M.Vab  = zeros(n,K);

M.Ea  = rand(n,K)*.25 + 1;      % make these positive
z = find(M.Ea<0);  M.Ea(z) = -M.Ea(z);
M.Eb  = randn(n,K)*.5;

M.Va  = rand(n,K)*2;      % make these positive
M.Vb  = rand(n,K)*2;      % make these positive
M.Vab = min(M.Va,M.Vb)./2;  % positive definite
M.Pik = exprnd(.5,n,K);
M.Pik = M.Pik ./ (sum(M.Pik,2)*ones(1,K));

% set the default maximum translation interval
if (isempty(M.Options.Interval))
  min_start = min(x(Seq(1:end-1)));
  max_end = max(x(Seq(2:end)-1));
  M.Options.Interval(2) = min_start-M.knots(1);
  M.Options.Interval(1) = max_end-M.knots(end);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SelectKnots
%%
function knots = SelectKnots(M,x,Y,Seq)
Start = min(x(Seq(1:end-1)));
End = max(x(Seq(2:end)-1));
Start = Start - M.Options.Interval(2);  % adjust for translation interval
End   = End   - M.Options.Interval(1);
len = max(diff(Seq)) + (M.Options.Interval(2)-M.Options.Interval(1));

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
Ops = SetFieldDef(Ops,'IterLimit',25);
Ops = SetFieldDef(Ops,'MaxDec',8);
Ops = SetFieldDef(Ops,'NumDec',0);
Ops = SetFieldDef(Ops,'NumEMStarts',1);
Ops = SetFieldDef(Ops,'MsgPrefix','');
Ops = SetFieldDef(Ops,'ShowGraphics',0);
Ops = SetFieldDef(Ops,'MinLen',[]);
Ops = SetFieldDef(Ops,'Interval',[]); % see InitE for default value
Ops = SetFieldDef(Ops,'Sigma',[]);
Ops.Sigma = SetFieldDef(Ops.Sigma,'Share',0);
Ops = SetFieldDef(Ops,'PropStart',0.8);
Ops = SetFieldDef(Ops,'InitializeWithExamples',0);
%
Ops.SearchOps = optimset('fminsearch');
Ops.SearchOps.Display = 'off';
Ops.SearchOps.MaxFunEvals = 200;
Ops.SearchOps.MaxIter = 200;

% Preprocessing
function Ops = showmodel(M,trajs)
M.Mu = permute(M.Mu, [1 3 2]);  % make p-K-D
[trash, M.C] = max(M.Pik,[],2);
Ops = viewmodel(M,trajs,M.Options);

function M = permuteModel(M)
M.Mu  = permute(M.Mu,[1 3 2]);

