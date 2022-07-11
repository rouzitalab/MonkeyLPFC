function M = lrm_ab(trajs,K,order,Ops)
%LRM_AB  LRM with transformation y=[ax+b]
%
%   Model = LRM_AB(Trajs,K,order,[Options])
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find
%    - order : order of polynomial regression
%
%   DefOps = LRM_AB('options) returns the default options.

% Scott Gaffney   1 Aug 2003
% Department of Information and Computer Science
% University of California, Irvine
%
% Changes
% -----------------

PROGNAME = 'lrm_ab';
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
[Y,x,Seq] = trajs2seq(trajs,M.zero,M.Options.MinLen);

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
M = InitE(M,Y,Seq);
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
  delete(Ops.MsgHnd);
end


%***************************************************************************
%   End Main Function
%***************************************************************************




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% postval 
%%
%%
function val = postval(pt,x,Y,Mu,order,r,s,sigma)
a = pt(1);  b = pt(2);
Xhat = regmat(a*x-b,order);
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
    n_i    = length(indx);
    
    % find the mode of the posterior and set Ea and Eb to this point
    a0 = M.Ea(i,k) - (1-M.Options.PropStart)*(M.Ea(i,k)-1);
    b0 = M.Eb(i,k) * M.Options.PropStart;
    pt0 = [a0 b0]';
    maxpt = ...
      fminsearch(fun,pt0,SearchOps,x(indx),Y(indx,:),Mu,P-1,r,s,sigma);
    a = maxpt(1);  b = maxpt(2);
    if (a==0), fprintf('a was zero.\n'); a=M.Ea(i,k); end
    M.Ea(i,k) = a;  M.Eb(i,k) = b;

    % Now we will calc the inverse Information at Ea and Eb.
    % First we make Xhat at (Ea,Eb) and then set up the gradients...
    Xhat = regmat(a*x(indx)-b,P-1);
    Dx = ones(n_i,1)*(0:P-1);
    Dx(:,3:end) = Dx(:,3:end) .* Xhat(:,2:end-1);
    D2x = ones(n_i,1)*[0 0 2:P-1];
    D2x(:,3:end) = D2x(:,3:end) .* Dx(:,2:end-1);
    
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
      Xhat = regmat(a(j)*x-b(j),M.order);
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
    fprintf(['lrm_ta_sh: Integration failed, using realmin*1e100 ',...
      'instead.\n']);
    zero = find(sum(M.Pik,2)==0);
    M.Pik(zero,:) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
    break;
  else
    fprintf(['lrm_ta_sh: Zero membership detected, trying ', ...
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
P = M.order+1;
K = M.K;
if (P==2)
  M = ExactMstep(M,x,Y,Seq);
  return;
end
[N,D] = size(Y);
n = length(Seq)-1;
lens = diff(Seq)';
MaxPD = 5;

% set up variance matrix indicies
rows  = 2;         diagi = [1:rows+1:rows^2];

%% Alpha
Wk = sum(M.Pik)';
ff = find(Wk==0);
if (~isempty(ff))
  fprintf('Some Wk are zero\n');
  Wk(ff) = 0.01*sum(Wk);
end
M.Alpha = Wk./n;

% initialization for Vxx and Y'*Vx
p = P-1;
Xp = cumprod([ones(N,1) x*ones(1,2*p)],2); % x^i => Xp(:,i+1)
Vxx = zeros(P,P);
YVx = zeros(P,D);
ga = zeros(n,2*p+1);  % g_i(a) => ga(:,i+1)
gb = zeros(n,2*p+1);  % g_i(b) => gb(:,i+1)
neg = -cumprod(-ones(1,2*p+1));  % neg->(0:2p) => neg(1:2p+1)

% initialization for covariance calculation
zz  = zeros(n,1);
m = 100;  % number of samples
a   = zeros(n,m,2*p);
b   = zeros(n,m,2*p);

% start the M-loop
for k=1:K
  Piik = copy(M.Pik(:,k),lens);

  %% s
  M.S(k) = M.Pik(:,k)'*(M.Eb(:,k).^2 + M.Vb(:,k));
  M.S(k) = M.S(k)./Wk(k);

  %% r
  M.R(k) = M.Pik(:,k)'*((M.Ea(:,k)-1).^2 + M.Va(:,k));
  M.R(k) = M.R(k)./Wk(k);


  % make power matrices for Ea,Eb,Va, and Vb
  Eap = cumprod([ones(n,1) M.Ea(:,k)*ones(1,2*p)],2); % a^i => Eap(:,i+1)
  Ebp = cumprod([ones(n,1) M.Eb(:,k)*ones(1,2*p)],2); % b^i => Ebp(:,i+1)
  Vap  = cumprod([M.Va(:,k)*ones(1,p)],2);   % Va^i => Vap(:,i)
  Vbp  = cumprod([M.Vb(:,k)*ones(1,p)],2);   % Vb^i => Vbp(:,i)
  
  % set up sampling to estimate cov(a^(p-i),b^i)
  for i=1:n
    Mn = [M.Ea(i,k); M.Eb(i,k)];
    V = [M.Va(i,k) M.Vab(i,k); M.Vab(i,k) M.Vb(i,k)];
    for j=1:MaxPD
      rnd = mvnrnd(Mn,V,m);
      if (isnan(rnd))
        if (j==MaxPD)
          fprintf('LRM_TA_SH: giving up, posterior covariance is not PD.\n');
          rnd = zeros(m,2); 
          break; 
        end
        V(diagi) = V(diagi)*1.5;   % bias trick, but alleviates probs
      else
        break;
      end
    end
    a(i,:,:) = cumprod(ones(2*p,1)*rnd(:,1)')';
    b(i,:,:) = cumprod(ones(2*p,1)*rnd(:,2)')';
  end

  % now compute YVx and Vxx
  for j=2:2*p
    fj = floor(j/2);
    MjCombo = ones(n,1)*(M.mj(1:fj).*M.Combo(j+1,2*(1:fj)+1));
    NegCombo = ones(N,1)*(neg(1:j+1).*M.Combo(j+1,1:j+1));

    % calculate gammas
    ga(:,j+1) = sum(MjCombo .* Vap(:,1:fj) .* Eap(:,j-1:-2:1),2);
    gb(:,j+1) = sum(MjCombo .* Vbp(:,1:fj) .* Ebp(:,j-1:-2:1),2);

    % calculate covariances
    if (j==2),  Vapbp = [zz M.Vab(:,k) zz];
    else,       Vapbp = [zz covv(a(:,:,j-1:-1:1),b(:,:,1:j-1)) zz];  end

    % calculate Gammas
    G_ab = Eap(:,j+1:-1:1) .*gb(:,1:j+1)    + ...
           Ebp(:,1:j+1)    .*ga(:,j+1:-1:1) + ...
           ga(:,j+1:-1:1)  .*gb(:,1:j+1)    + Vapbp;

    % calculate Deltas
    Vx = Piik.*sum(NegCombo.*copy(G_ab,lens).*Xp(:,j+1:-1:1),2);
    if (j<=p)
      YVx(j+1,:) = Vx'*Y;
    end

    % place the results into the matrix form of Vxx
    p1=1; p2=j+1;
    dj = p2-P;  % make sure Vxx is big enough for the whole diagonal
    if (dj>0), p1=p1+dj;  p2=p2-dj; end % if not, then trim it down
    Vxx(P.*([p2:-1:p1]-1) + [p1:p2]) = sum(Vx); % copy correct diagonal
  end


  % now make Xhat and finish up the M-Step
  Xhat = regmat(copy(M.Ea(:,k),lens).*x- copy(M.Eb(:,k),lens),P-1);
  PikXhat = Piik*ones(1,P).*Xhat;   

  %% Mu
  M.Mu(:,:,k) = (PikXhat'*Xhat + Vxx) \ (PikXhat'*Y + YVx);

  % Sigma
  M.Sigma(k,:) = (Piik'*((Y-Xhat*M.Mu(:,:,k)).^2) ...
    -2*sum(YVx.*M.Mu(:,:,k)) + sum((M.Mu(:,:,k)'*Vxx)'.*M.Mu(:,:,k))) ...
    ./ sum(M.Pik(:,k).*lens);
end
if (M.Options.Sigma.Share)
  M.Sigma = (sum(M.Sigma,2)./D)*ones(1,D);
end

% check for strange circumstances
fr = find(M.R<=0);  fs = find(M.S<=0);
fsig = find(M.Sigma<=0);
if (~isempty([fr(:); fs(:); fsig(:)]))
  fprintf('Some variances are zero\n');
  M.R(fr) = M.Options.minvar;  M.S(fs) = M.Options.minvar;
  M.Sigma(fsig) = M.Options.minvar;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ExactMstep
%%
function M = ExactMstep(M,x,Y,Seq)
[N,D] = size(Y);
P = M.order+1;
K = M.K;
n = length(Seq)-1;
lens = diff(Seq)';
Wk = sum(M.Pik)';
ff = find(Wk==0);
if (~isempty(ff))
  fprintf('Some Wk are zero\n');
  Wk(ff) = 0.01*sum(Wk);
end

%% Alpha
M.Alpha = Wk./n;

% hyper prior stuff
nu_r = 10;
r0   = (nu_r-2)/nu_r *.1;   % E[M.R] = .1
nu_s = 10;
s0   = (nu_s-2)/nu_s *.3;   % E[M.S] = .3
nu_r=2; nu_s=2; r0=0; s0=0;   % erase prior w/ this line

for k=1:K
  Piik = copy(M.Pik(:,k),lens);
  
  %% s
  M.S(k) = M.Pik(:,k)'*(M.Eb(:,k).^2 + M.Vb(:,k));
  M.S(k) = (nu_s*s0 + M.S(k))./ (Wk(k) + nu_s-2);
  
  %% r
  M.R(k) = M.Pik(:,k)'*((M.Ea(:,k)-1).^2 + M.Va(:,k));
  M.R(k) = (nu_r*r0 + M.R(k))./ (Wk(k) + nu_r-2);
  
  % prepare for Mu and Sigma
  trLam = sum(Piik.*x.^2.*copy(M.Va(:,k),lens)) + ...
    sum(lens.*M.Pik(:,k).*M.Vb(:,k)) - ...
    2*sum(Piik.*x.*copy(M.Vab(:,k),lens));
  Xhat = regmat(copy(M.Ea(:,k),lens).*x-copy(M.Eb(:,k),lens),P-1);
  PikXhat = Piik*ones(1,P).*Xhat;   

  %% Mu
  M.Mu(:,:,k) = (PikXhat'*Xhat + [0 0; 0 trLam]) \ (PikXhat'*Y);
    
  %% Sigma
  M.Sigma(k,:) = ( sum(Piik'*((Y-Xhat*M.Mu(:,:,k)).^2)) + ...
    M.Mu(2,:,k)^2*trLam ) ./ sum(lens.*M.Pik(:,k));
end
if (M.Options.Sigma.Share)
  M.Sigma = (sum(M.Sigma,2)./D)*ones(1,D);
end

% check for strange circumstances
fr = find(M.R<=0);  fs = find(M.S<=0);
fsig = find(M.Sigma<=0);
if (~isempty([fr(:); fs(:); fsig(:)]))
  fprintf('Some variances are zero\n');
  M.R(fr) = M.Options.minvar;  M.S(fs) = M.Options.minvar;
  M.Sigma(fsig) = M.Options.minvar;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% InitM
%%
function M = InitM(M,x,Y,Seq)
M = Mstep(M,x,Y,Seq);

if (M.Options.InitializeWithExamples)
[P,D,K] = size(M.Mu);
X = regmat(x,P-1);
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
function M = InitE(M,Y,Seq)
[N,D] = size(Y);
K = M.K;
P = M.order+1;
n = length(Seq)-1;

% Set up combination matrix and variance coefficient vector used in the M-step
maxp = 2*(P-1);
M.Combo = ones(maxp+1,maxp+1);
for p=1:maxp
  for q=1:p-1
    if (q>floor(p/2))
      M.Combo(p+1,q+1) = M.Combo(p+1,p-q+1);
    else
      M.Combo(p+1,q+1) = prod(p-q+1:p)/prod(1:q);
    end
  end
end
M.mj = cumprod(1:2:maxp);  % product of odd numbers


% E-step vars
% M.Ea   = ones(n,K);
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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%
% Options Handling

function Ops = DefaultOptions(Ops);
Ops = SetFieldDef(Ops,'zero','nozero');
Ops = SetFieldDef(Ops,'stopval',1e-7);
Ops = SetFieldDef(Ops,'minvar',1e-5);
Ops = SetFieldDef(Ops,'MaxDec',4);
Ops = SetFieldDef(Ops,'NumDec',0);
Ops = SetFieldDef(Ops,'IterLimit',25);
Ops = SetFieldDef(Ops,'NumEMStarts',1);
Ops = SetFieldDef(Ops,'MsgPrefix','');
Ops = SetFieldDef(Ops,'ShowGraphics',0);
Ops = SetFieldDef(Ops,'MsgHnd',[]);
Ops = SetFieldDef(Ops,'MinLen',[]);
Ops = SetFieldDef(Ops,'Sigma',[]);
Ops.Sigma = SetFieldDef(Ops.Sigma,'Share',0);
Ops = SetFieldDef(Ops,'PropStart',0.8);
Ops = SetFieldDef(Ops,'InitializeWithExamples',1);
%
Ops.SearchOps = optimset('fminsearch');
Ops.SearchOps.Display = 'off';
Ops.SearchOps.MaxFunEvals = 200;
Ops.SearchOps.MaxIter = 200;


function c = covv(x,y)
[m,n,p] = size(x);
c = permute((sum(x.*y,2) - sum(x,2).*sum(y,2)./n)./(n),[1 3 2]);

function Ops = showmodel(M,trajs)
M.Mu = permute(M.Mu, [1 3 2]);  % make p-K-D
[trash, M.C] = max(M.Pik,[],2);
Ops = viewmodel(M,trajs,M.Options);

function M = permuteModel(M)
M.Mu  = permute(M.Mu,[1 3 2]);

