function M = srm_d_b(trajs,K,knots,order,Ops)
%SRM_D_B  SRM with transformation y=[x+b]+d
%
%   Model = SRM_D_B(Trajs,K,knots,order,[Options])
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find
%    - knots : vector of knot values (pass in [] if you want the algorithm
%              to select the knots "automatically"; Use Options.Kn if
%              you would like to specify the number of knots); note that
%              the resulting augmented knot sequence may have a longer
%              length (this is normal).
%    - order : order of polynomial regression
%
%   DefOps = SRM_D_B('options) returns the default options.
%
%   Options (some fields)
%    .Interval  : range of allowed translations (e.g., [-2 2]), you  must
%                  provide this if knots is empty when this functions is called.
%    .KnotMult : (default .25) multiplied by the complete length of data 
%                interval to get the number of knots (minimum of 1)

% Scott Gaffney   9 October 2003
% Department of Information and Computer Science
% University of California, Irvine
%
% Changes
% -----------------

PROGNAME = 'srm_d_b';
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
    errorbox([PROGNAME,': you must provide either ', ...
        '''knots'' or ''Options.Interval''.']);
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
  M = Estep(M,x,Y,Seq);
  [Lhood(NumIter), M.Pik] = CalcLike(M.Pik,PROGNAME);
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
M.NumIndParams = (K-1) + K*P*D + 2*K*D + K;  % alpha, mu, sigma, t, s
if (M.Options.Sigma.Share==1)
  M.NumIndParams = M.NumIndParams - K*(D-1);
end
if (M.Options.Tf.Share==1)
  M.NumIndParams = M.NumIndParams - K*(D-1);
end

if (CreatedMsgBar)
  delete(Ops.MsgHnd);
end


%***************************************************************************
%   End Main Function
%***************************************************************************




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% postval 
%%
%%
function val = postval(pt,x,Y,Mu,knots,order,s,t,sigma)
b = pt(1); f = pt(2:end)';
n = length(x);
Xhat = bsplinebasis(knots,order,x-b);
val = sum(sum((Y-Xhat*Mu-f(ones(n,1),:)).^2)./sigma) + sum(f.^2./t) + b^2/s;

function val = b_y(pt,x,Y,Mu,knots,order,s,t,sigma)
b = pt(1);
ni = length(x);
D = length(t);
val = 0;
Xhat = bsplinebasis(knots,order,x-b);
YxMu = Y-Xhat*Mu;
for d=1:D
  iV = eye(ni)./sigma(d) - 1/(ni*sigma(d) + sigma(d)^2/t(d));
  val = val + YxMu(:,d)'*iV*YxMu(:,d);
end
val = val + b^2/s;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estep
%%
function M = Estep(M,x,Y,Seq)
[P,D,K] = size(M.Mu);
n = length(Seq)-1;
fun = @postval;
bfun = @b_y;
SearchOps = M.Options.SearchOps;

% find posterior mode for Eb and Ef
for k=1:K
  t        = M.T(k,:);
  s        = M.S(k);
  Mu       = M.Mu(:,:,k);
  sigma    = M.Sigma(k,:);
  for i=1:n
    indx   = Seq(i):Seq(i+1)-1;
    ni    = length(indx);
    b0 = M.Eb(i,k) * M.Options.PropStart;
    pt0 = [b0];
    b = fminsearch(bfun,pt0,SearchOps,x(indx),Y(indx,:),Mu, ...
      M.knots,M.order,s,t,sigma);
    Xhat = bsplinebasis(M.knots,M.order,x(indx)-b);
    f  = t/(ni*t+sigma)*sum(Y(indx,:)-Xhat*Mu);
    M.Eb(i,k) = b;   M.Ef(i,:,k) = f;
    
    % Now we will calc the inverse Information at Eb and Ef.
    % First we calculate the derivatives.
    Dx  = bsplinebasis(M.knots,M.order,x(indx)-b,1);
    D2x = bsplinebasis(M.knots,M.order,x(indx)-b,2);
    
    % ...then make repeat calculations
    DxMu = Dx*Mu;  % ni x D
    D2xMu = D2x*Mu;  % ni x D
    YxMu = Y(indx,:)-Xhat*Mu-f(ones(ni,1),:);  % ni x D
    
    % ...and now we can calculate the information.
    Ib  = sum(sum(YxMu.*D2xMu)./sigma) - sum(sum(DxMu.^2)./sigma) - 1./s;
    If  = -ni./sigma - 1./t;
    Ibf = sum(DxMu)./sigma;
    
    % the variance is the inverse of the negative information
    V = diag([Ib If]);
    V(1,2:end) = Ibf;  V(2:end,1) = Ibf';
    V = inv(-V);
    var = diag(V);
    M.Vb(i,k)    = var(1);
    M.Vf(i,:,k)  = var(2:end);
    M.Vbf(i,:,k) = V(1,2:end);
    %M.Vbf(i,:,k) = 0;
  end
end

M = CalcPik(M,x,Y,Seq);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcPik 
%%
function M = CalcPik(M,x,Y,Seq)
% Numerical integration
NumSamps = 70;
MaxTries = 5;
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
  M.Pik = (Pik./TotalSamps) .* (ones(n,1)*M.Alpha'); % we keep Pik for next try!
  if (all(sum(M.Pik,2))), break; end

  % we have detected some zeros, try again?
  if (tries==MaxTries)
    fprintf('Integration failed, using realmin*1e100 instead.\n');
    zero = find(sum(M.Pik,2)==0);
    M.Pik(zero,:) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
    break;
  else
    fprintf('Zero membership detected, trying integration again: %d\n',tries);
    tries = tries+1;
    %NumSamps = round(NumSamps*1.25);
    S = 1.25*S;  % biased, but gets over some tricky integrations
    %T = 1.25*T;
  end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcLike
%%
function [Lhood, Pik] = CalcLike(Pik,PROGNAME)
[n,K] = size(Pik);
s = sum(Pik,2);
Lhood = sum(log(s));
Pik = Pik ./ (s*ones(1,K));  % normalize the memberships while were at it




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% StoppingCondition
%%
function [dostop,Ops] = StoppingCondition(Lhood,NumIter,Ops)
dostop=0;

if (NumIter >= Ops.IterLimit)
  dostop = 1;
  return;
end
% return;  % for alignment

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
function M = Mstep(M,x,Y,Seq)
K = M.K;
[N,D] = size(Y);
n = length(Seq)-1;
P = length(M.knots)-M.order;
lens = diff(Seq)';
Wk = sum(M.Pik)';

%% Alpha
M.Alpha = Wk./n;

% start the M-loop
for k=1:K
  %% s
  M.S(k) = M.Pik(:,k)'*(M.Eb(:,k).^2 + M.Vb(:,k));
  M.S(k) = M.S(k)./Wk(k);

  %% r
  M.T(k,:) = M.Pik(:,k)'*(M.Ef(:,:,k).^2 + M.Vf(:,:,k));
  M.T(k,:) = M.T(k,:)./Wk(k);
  
  %make Xhat and finish up the M-Step
  Xhat = bsplinebasis(M.knots,M.order,x-copy(M.Eb(:,k),lens));
  Piik = copy(M.Pik(:,k),lens);
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
  M.valid(:,k) = [valid(1);valid(end)];
  PikXhat(:,invalid) = [];
  
  %% Mu
  Ef = copy(M.Ef(:,:,k),lens);
  M.Mu(:,:,k) = (PikXhat'*Xhat) \ (PikXhat'*(Y-Ef));
  
  % fill-in invalid coefficients w/ closest valid ones
  f = find(invalid<valid(1));
  M.Mu(invalid(f),:,k) = ones(length(f),1) * M.Mu(valid(1),:,k);
  f = find(invalid>valid(end));
  M.Mu(invalid(f),:,k) = ones(length(f),1) * M.Mu(valid(end),:,k);

  % Sigma
  M.Sigma(k,:) = ( Piik'*((Y-Xhat*M.Mu(:,:,k)-Ef).^2) ...
    + (M.Pik(:,k).*lens)'*M.Vf(:,:,k) )./ sum(M.Pik(:,k).*lens);
end
if (M.Options.Sigma.Share)
  M.Sigma = (sum(M.Sigma,2)./D)*ones(1,D);
end
if (M.Options.Tf.Share)
  M.T = (sum(M.T,2)./D)*ones(1,D);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% InitM
%%
function M = InitM(M,x,Y,Seq)
M = Mstep(M,x,Y,Seq);  % default init

if (M.Options.InitializeWithExamples)
% this is key for this method
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
M.Eb   = zeros(n,K);
M.Ef   = zeros(n,D,K);
M.Vb   = zeros(n,K);
M.Vf   = zeros(n,D,K);
M.Vbf  = zeros(n,D,K);

M.Eb  = randn(n,K)*.5;
M.Ef  = randn(n,D,K)*.5;
M.Vb  = rand(n,K)*2;        % make these positive
M.Vf  = rand(n,D,K)*2;      % make these positive
tempVb = permute(M.Vb(:,:,ones(1,D)),[1 3 2]);
M.Vbf = min(tempVb,M.Vf)./2;  % positive definite
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
  Kn = max(ceil(len/3),1);
end
knots = linspace(Start,End,Kn);  % uniform spacing
knots = addendpts(knots,M.order);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%

function Ops = DefaultOptions(Ops);
Ops = SetFieldDef(Ops,'zero','nozero');
Ops = SetFieldDef(Ops,'stopval',1e-7);
Ops = SetFieldDef(Ops,'minvar',1e-5);
Ops = SetFieldDef(Ops,'IterLimit',15);
Ops = SetFieldDef(Ops,'MaxDec',4);
Ops = SetFieldDef(Ops,'NumDec',0);
Ops = SetFieldDef(Ops,'NumEMStarts',1);
Ops = SetFieldDef(Ops,'MsgPrefix','');
Ops = SetFieldDef(Ops,'ShowGraphics',0);
Ops = SetFieldDef(Ops,'MinLen',[]);
Ops = SetFieldDef(Ops,'Interval',[]); % see InitE for default value
Ops = SetFieldDef(Ops,'PropStart',0.8);
Ops = SetFieldDef(Ops,'Sigma',[]);
Ops.Sigma = SetFieldDef(Ops.Sigma,'Share',0);
Ops = SetFieldDef(Ops,'Tf',[]);
Ops.Tf = SetFieldDef(Ops.Tf,'Share',0);
Ops = SetFieldDef(Ops,'InitializeWithExamples',0);
%
Ops.SearchOps = optimset('fminsearch');
Ops.SearchOps.Display = 'off';
Ops.SearchOps.MaxFunEvals = 200;
Ops.SearchOps.MaxIter = 200;

function Ops = showmodel(M,trajs)
M = permuteModel(M);
[trash, M.C] = max(M.Pik,[],2);
Ops = viewmodel(M,trajs,M.Options);

function M = permuteModel(M)
M.Mu  = permute(M.Mu,[1 3 2]);
M.Ef  = permute(M.Ef,[1 3 2]);

