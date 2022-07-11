function [M,Ops] = lrm_b(trajs,K,order,Ops)
%LRM_B  LRM with transformation y=[x+b]
%
%   Model = LRM_B(Trajs,K,order,[Options])
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find
%    - order : order of polynomial regression
%
%   DefOps = LRM_B('options) returns the default options.

% Scott Gaffney   8 May 2003
% Department of Information and Computer Science
% University of California, Irvine
%
% Changes
% -----------------

PROGNAME = 'lrm_b';
METHOD   = PROGNAME;
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
  [Lhood(NumIter),M] = CalcLike(M,N);
  [DoStop,Ops] = StoppingCondition(Lhood,NumIter,Ops);
  if (DoStop), break; end

  %%% M-Step
  M = Mstep(M,x,Y,Seq);
end
%%%%%%%%%%%%%%%%%%% E-M Algorithm

M.NumPoints = prod(size(Y));
M = permuteModel(M); 
M.Lhood = Lhood(1:NumIter);
[trash, M.C] = max(M.Pik,[],2);
M.TrainLhood = M.Lhood(end);
M.TrainLhood_ppt = M.Lhood(end)./M.NumPoints;

% Calculate number of independent parameters
[P,K,D] = size(M.Mu);
M.NumIndParams = (K-1) + K*P*D + K*D + K;  % alpha, mu, sigma, s
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
function val = postval(b,x,Y,Mu,order,sigma,s)
Xhat = regmat(x-b,order);
val = sum(sum((Y-Xhat*Mu).^2)./sigma) + (b.^2)./s;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estep 
%%
function M = Estep(M,x,Y,Seq)
[P,D,K] = size(M.Mu);
[N,D] = size(Y);
n = length(Seq)-1;
lens = diff(Seq);
fun = @postval;
SearchOps = M.Options.SearchOps;

%%%% Calculate posterior mode for p(delta|y) and evaluate 
%%%%   V[delta|y] at that value
if (P>2)
  for k=1:K
    Mu    = M.Mu(:,:,k);
    sigma = M.Sigma(k,:);
    s     = M.S(k);
    
    for i=1:n
      indx = Seq(i):Seq(i+1)-1;
      ni = length(indx);
      b0 = M.Eb(i,k) * M.Options.PropStart;
      pt0 = [b0];
      b =  fminsearch(fun,pt0,SearchOps,x(indx),Y(indx,:),Mu,P-1,sigma,s);
      M.Eb(i,k) = b;
      
      % Now we will calc the inverse Information at Eb.
      % First we make Xhat at Eb and then set up the derivatives
      Xhat = regmat(x(indx)-M.Eb(i,k),P-1);
      Dx = ones(ni,1)*(0:P-1);
      Dx(:,3:end) = Dx(:,3:end) .* Xhat(:,2:end-1);
      D2x = ones(ni,1)*[0 0 2:P-1];
      D2x(:,3:end) = D2x(:,3:end) .* Dx(:,2:end-1);
      
      % calculate inverse information
      DxMu = Dx*Mu;
      YxMu = Y(indx,:)-Xhat*Mu;
      Ib = sum( (sum(YxMu.*(D2x*Mu)) - sum(DxMu.^2))./sigma ) -1/s;
      M.Vb(i,k) = -1/Ib;
    end
  end
  
% Exact case for order 1 regression models
else
  X = [ones(N,1) x];
  for k=1:K
    M.Vb(:,k) = 1./(sum((M.Mu(2,:,k).^2)./M.Sigma(k,:)).*lens + 1/M.S(k));
    YXMuk = Y-X*M.Mu(:,:,k);
    for i=1:n
      indx = Seq(i):Seq(i+1)-1;
      M.Eb(i,k) = ...
        M.Vb(i,k)*sum((-M.Mu(2,:,k)./M.Sigma(k,:)).*sum(YXMuk(indx,:)));
    end
  end
end


if (P>2)
  M = CalcPik(M,x,Y,Seq);
else
  X = [ones(N,1) x];
  Pikd = zeros(n,K,D);
  for k=1:K
    br = M.Mu(2,:,k).^2*M.S(k);
    YXMuk = Y-X*M.Mu(:,:,k);
    for i=1:n
      for d=1:D
        indx = Seq(i):Seq(i+1)-1;
        ni = length(indx);
        Vi = br(d)*ones(ni) + M.Sigma(k,d)*eye(ni);
        Pikd(i,k,d) = M.Alpha(k)*mvnormpdf(YXMuk(indx,d)',zeros(ni,1),Vi);
      end
    end
  end
  M.Pik = prod(Pikd,3);
end


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
P = M.order+1;
mlen = max(diff(Seq));
Piid = zeros(N,D);
Piimk = zeros(N,NumSamps,K);
S = M.S;
TotalSamps = 0;
tries = 1;

while (1)
  TotalSamps = TotalSamps + NumSamps;
  M.Pik(:) = 0;

  % calculate the density at sampled points
  for k=1:K
    b = randn(NumSamps,1).*sqrt(S(k));
    for j=1:NumSamps
      Xhat = regmat(x-b(j),P-1);
      for d=1:D
        Piid(:,d) = normpdf(Y(:,d),Xhat*M.Mu(:,d,k),M.Sigma(k,d));
      end
      Piimk(:,j,k) = prod(Piid,2);  % get rid of D here (mult it out)
    end
  end

  % now scale the data to avoid underflow with long curves
  % and sum across the sample integration points
  M.scale = mean(mean(mean(Piimk)));
  Piimk_scl = Piimk./M.scale;  % we don't scale across D; we got rid of D above.
  for k=1:K
    for j=1:TotalSamps
      M.Pik(:,k) = M.Pik(:,k) + sprod(Piimk_scl(:,j,k),Seq,mlen);
    end
  end
  clear Piimk_scl;
  M.Pik = (M.Pik./TotalSamps) .* (ones(n,1)*M.Alpha');
  if (all(sum(M.Pik,2))), break; end  % check for stable integration

  % we have detected some zeros, try again?
  if (tries==MaxTries)
    fprintf('lrm_tt_sh: Integration failed, using realmin*1e100 instead.\n');
    zero = find(sum(M.Pik,2)==0);
    M.Pik(zero,:) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
    break;
  else
    fprintf(['lrm_tt_sh: Zero membership detected, trying integration ',...
      'again: %d\n'],tries);
    tries = tries+1;
    S = 1.25*S;  % biased, but gets over some tricky integrations
    Piimk = [zeros(N,NumSamps,K) Piimk]; % save current values
  end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcLike
%%
function [Lhood,M] = CalcLike(M,N)
K = M.K;
P = size(M.Mu,1);
s = sum(M.Pik,2);
Lhood = sum(log(s));
if (P>2),  Lhood = Lhood + N*log(M.scale); end
M.Pik = M.Pik ./ (s*ones(1,K));  % normalize the memberships




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

% Alpha
M.Alpha = Wk./n;

% initialization for Vxx and Y'*Vx
p = P-1;
Xp = cumprod([ones(N,1) x*ones(1,2*p)],2); % x^i => Xp(:,i+1)
Vxx = zeros(P,P);
YVx = zeros(P,D);
gb = zeros(n,2*p+1);  % g_i(b) => gb(:,i+1)
neg = -cumprod(-ones(1,2*p+1));  % neg->(0:2p) => neg(1:2p+1)

% hyper prior stuff
Es = 4;  % E[M.S] = Es
nu_s = 10;
s0   = (nu_s-2)/nu_s *Es;   
nu_s=2; s0=0;   % erase prior w/ this line


% start the M-loop
for k=1:K
  Piik = copy(M.Pik(:,k),lens);

  % s
  M.S(k) = M.Pik(:,k)'*(M.Eb(:,k).^2 + M.Vb(:,k));
  M.S(k) = (nu_s*s0 + M.S(k))./ (Wk(k) + nu_s-2);
  
  % make power matrices for Eb and Vb
  Ebp = cumprod([ones(n,1) M.Eb(:,k)*ones(1,2*p)],2); % b^i => Ebp(:,i+1)
  Vbp  = cumprod([M.Vb(:,k)*ones(1,p)],2);   % Vb^i => Vbp(:,i)

  % now compute YVx and Vxx
  for j=1:2*p
    fj = floor(j/2);
    MjCombo = ones(n,1)*(M.mj(1:fj).*M.Combo(j+1,2*(1:fj)+1));
    NegCombo = ones(N,1)*(neg(1:j+1).*M.Combo(j+1,1:j+1));

    % calculate gammas and Gammas
    gb(:,j+1) = sum(MjCombo.* Vbp(:,1:fj) .* Ebp(:,j-1:-2:1),2);
    Gb = gb(:,1:j+1);
    
    % calculate Deltas
    Vx = Piik.*sum(NegCombo.*copy(Gb,lens).*Xp(:,j+1:-1:1),2);
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
  Xhat = regmat(x-copy(M.Eb(:,k),lens),P-1);
  PikXhat = Piik*ones(1,P).*Xhat;   

  % Mu
  M.Mu(:,:,k) = (PikXhat'*Xhat + Vxx) \ (PikXhat'*Y + YVx);

  M.Sigma(k,:) = ( Piik'*((Y-Xhat*M.Mu(:,:,k)).^2) ...
    -2*sum(YVx.*M.Mu(:,:,k)) + sum((M.Mu(:,:,k)'*Vxx)'.*M.Mu(:,:,k)) ) ...
    ./ sum(M.Pik(:,k).*lens);
end
if (M.Options.Sigma.Share)
  M.Sigma = (sum(M.Sigma,2)./D)*ones(1,D);
end

% check for strange circumstances
fs = find(M.S<=0);
fsig = find(M.Sigma<=0);
if (~isempty([fs(:); fsig(:)]))
  fprintf('Some variances are zero\n');
  M.S(fs) = M.Options.minvar;
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
M.Eb   = zeros(n,K);
M.Vb   = zeros(n,K);

M.Eb  = randn(n,K)*.5;
M.Vb  = rand(n,K)*2;
M.Pik = exprnd(.5,n,K);
M.Pik = M.Pik ./ (sum(M.Pik,2)*ones(1,K));





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
%%

function Ops = DefaultOptions(Ops);
Ops = SetFieldDef(Ops,'zero','nozero');
Ops = SetFieldDef(Ops,'stopval',1e-5);
Ops = SetFieldDef(Ops,'minvar',1e-5);
Ops = SetFieldDef(Ops,'MaxDec',4);
Ops = SetFieldDef(Ops,'NumDec',0);
Ops = SetFieldDef(Ops,'IterLimit',50);
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

function Ops = showmodel(M,trajs)
M = permuteModel(M);
[trash, M.C] = max(M.Pik,[],2);
Ops = viewmodel(M,trajs,M.Options);

function M = permuteModel(M)
M.Mu  = permute(M.Mu,[1 3 2]);

