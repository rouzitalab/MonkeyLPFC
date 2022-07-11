function M = lrm_cd_b(trajs,K,order,Ops)
%LRM_CD_B  LRM with transformation y=c[x+b]+d
%
%   Model = LRM_CD_B(Trajs,K,order,[Options])
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find
%    - order : order of polynomial regression
%
%   DefOps = LRM_CD_B('options) returns the default options.

% Scott Gaffney   1 Aug 2003
% Department of Information and Computer Science
% University of California, Irvine
%
% Changes
% -----------------

PROGNAME = 'lrm_cd_b';
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
M.Options.ScaleInSpace=1;  % this is used for plotting, e.g. see Plot_Regmix()

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
  [Lhood(NumIter),M] = CalcLike(M);
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
M.NumIndParams = (K-1) + K*P*D + 3*K*D + K;  % alpha, mu, sigma, u, t, s
if (M.Options.Sigma.Share==1)
  M.NumIndParams = M.NumIndParams - K*(D-1);
end
if (M.Options.Tf.Share==1)
  M.NumIndParams = M.NumIndParams - K*(D-1);
end
if (M.Options.U.Share==1)
  M.NumIndParams = M.NumIndParams - K*(D-1);
end

if (CreatedMsgBar)
  delete(Options.MsgHnd);
end


%***************************************************************************
%   End Main Function
%***************************************************************************




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% postval 
%%
%%
function val = b_y(pt,x,Y,Mu,order,s,u,t,sigma)
b = pt(1);
ni = length(x);
D = length(t);
val = 0;
Xhat = regmat(x-b,order);
XMu = Xhat*Mu;
YxMu = Y-XMu;
for d=1:D
  iR = eye(ni)/sigma(d) - XMu(:,d)*XMu(:,d)'/(sigma(d)^2/u(d) ...
    + sigma(d)*(XMu(:,d)'*XMu(:,d)));
  siR = sum(iR);
  iV = iR - sum(iR,2)*siR/(1/t(d) + sum(siR));
  val = val + YxMu(:,d)'*iV*YxMu(:,d);
end
val = val + b^2/s;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estep
%%
function M = Estep(M,x,Y,Seq)
[P,D,K] = size(M.Mu);
n = length(Seq)-1;
fun = @b_y;
pfun = @post;
SearchOps = M.Options.SearchOps;

% set up variance matrix indicies
rows  = 1+2*D;         diagi = [1:rows+1:rows^2];
rws_e = 2:2+D-1;       rws_f = 2+D:2+D+D-1;
bb_ij = diagi(1);
ee_ij = diagi(rws_e);  ff_ij = diagi(rws_f);
ef_ij = rws_e + rows*(rws_f-1);
be_ij = rws_e+rows;    bf_ij = rws_f+rows;

for k=1:K
  s        = M.S(k);
  u        = M.U(k,:);
  t        = M.T(k,:);
  Mu       = M.Mu(:,:,k);
  sigma    = M.Sigma(k,:);
  for i=1:n
    indx   = Seq(i):Seq(i+1)-1;
    ni    = length(indx);  

    b0 = M.Eb(i,k) * M.Options.PropStart;
    pt0 = [b0];
    M.Eb(i,k) = fminsearch(fun,pt0,SearchOps, ...
                           x(indx),Y(indx,:),Mu,P-1,s,u,t,sigma);

    % Find e and f at b_hat
    Xhat = regmat(x(indx)-M.Eb(i,k),P-1);
    XMu  = Xhat*Mu;
    I    = eye(ni);
    MuXXMu = sum(XMu.*XMu);
    for d=1:D
      siR = sum(I/sigma(d) - ...
        XMu(:,d)*XMu(:,d)'/(sigma(d)^2/u(d) + sigma(d)*MuXXMu(d)));
      iS = I/sigma(d) - 1/(ni*sigma(d) + sigma(d)^2/t(d));
      MuXiS = (XMu(:,d)'* iS);
      Ve  = u(d)/ (u(d)*MuXiS* XMu(:,d) + 1);
      Vf  = t(d)/ (t(d)*sum(siR) + 1);
      M.Ee(i,d,k)  = Ve*(MuXiS*Y(indx,d) +1/u(d));
      M.Ef(i,d,k)  = Vf*siR*(Y(indx,d)-XMu(:,d));
    end

    % Now we will calc the inverse Information at Eb and Ee,Ef.
    Dx = ones(ni,1)*(0:P-1);
    Dx(:,3:end) = Dx(:,3:end) .* Xhat(:,2:end-1);
    D2x = ones(ni,1)*[0 0 2:P-1];
    D2x(:,3:end) = D2x(:,3:end) .* Dx(:,2:end-1);
    
    % ...then make repeat calculations
    Mue = Mu.*(ones(P,1)*M.Ee(i,:,k));
    f = ones(ni,1)*M.Ef(i,:,k);
    Xmu     = Xhat*Mu; % ni x D
    Dxmu    = Dx*Mu;  % ni x D
    DxMu    = Dx*Mue;  % ni x D
    D2xMu   = D2x*Mue;  % ni x D
    YxMu = Y(indx,:) - Xhat*Mue - f;    % ni x D

    % ...and now we can calculate the information.
    Ib  = sum(sum(YxMu.*D2xMu)./sigma)   - sum(sum(DxMu.^2)./sigma)  - 1./s;
    Ibe = -(sum(YxMu.*Dxmu) - sum(Xmu.*Dxmu))./sigma;
    Ibf = sum(DxMu)./sigma;
    Ie  = -sum(Xmu.^2)./sigma - 1./u;
    If  = -ni./sigma - 1./t;
    Ief = -sum(Xmu)./sigma;
    
    % the variance is the inverse of the negative information
    V = diag([Ib Ie If]);
    V(1,rws_e) = Ibe;  V(rws_e,1) = Ibe';
    V(1,rws_f) = Ibf;  V(rws_f,1) = Ibf';
    V(ef_ij)   = Ief;  V(ef_ij-D*(rows-1)) = Ief;

    V = inv(-V);
    M.Vb(i,k)    = V(bb_ij);
    M.Ve(i,:,k)  = V(ee_ij);
    M.Vf(i,:,k)  = V(ff_ij);
    M.Vef(i,:,k) = V(ef_ij);
    M.Vbe(i,:,k) = V(be_ij);
    M.Vbf(i,:,k) = V(bf_ij);
  end
end
fb = find(M.Vb<=0);
fe = find(M.Ve<=0);    ff = find(M.Vf<=0);
if (~isempty([fb(:); fe(:); ff(:)]))
  fprintf('A negative variance was detected....\n');
  M.Vb(fb) = -M.Vb(fb);
  M.Ve(fe) = -M.Ve(fe);  M.Vf(ff) = -M.Vf(ff);
end

M = CalcPik(M,x,Y,Seq);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcPik 
%%
function M = CalcPik(M,x,Y,Seq)
% Numerical integration
NumSamps = 150;
MaxTries = 5;
[N,D] = size(Y);
n = length(Seq)-1;
P = M.order+1;
K = M.K;
Pid = zeros(n,D);
Pik = zeros(n,K);
S = M.S;  U = M.U;  T = M.T;
TotalSamps = 0;
tries = 1;

while (1)
  TotalSamps = TotalSamps + NumSamps;

  % calculate the density at sampled points
  for k=1:K
    b = randn(NumSamps,1).*sqrt(S(k));   % sample from N(0,s)
    for j=1:NumSamps
      Xhat = regmat(x-b(j),P-1);
      XMu = Xhat*M.Mu(:,:,k);
      for d=1:D
        sigma = M.Sigma(k,d);  u = U(k,d);  t = T(k,d);
        for i=1:n
          indx = Seq(i):Seq(i+1)-1;  ni = length(indx);
          iR = eye(ni)./sigma - XMu(indx,d)*XMu(indx,d)'/(sigma^2/u ...
            + sigma*(XMu(indx,d)'*XMu(indx,d)));  siR = sum(iR);
          iV = iR - sum(iR,2)*siR/(1/t + sum(siR));
          Pid(i,d) = mvnormpdf_inv(Y(indx,d)',XMu(indx,d),iV);
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
    %S = 1.25*S;  % biased, but gets over some tricky integrations
    %R = 1.25*T;
  end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcLike
%%
function [Lhood, M] = CalcLike(M)
[n,K] = size(M.Pik);
s = sum(M.Pik,2);
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
function M = Mstep(M,x,Y,Seq)
P = M.order+1;
K = M.K;
[N,D] = size(Y);
n = length(Seq)-1;
lens = diff(Seq)';
Wk = sum(M.Pik)';
ff = find(Wk==0);
if (~isempty(ff))
  fprintf('Some Wk are zero\n');
  Wk(ff) = 0.01*sum(Wk);
end
MaxPD = 6;

% set up variance matrix indicies
rows  = 1+2*D;         diagi = [1:rows+1:rows^2];
rws_e = 2:2+D-1;       rws_f = 2+D:2+D+D-1;
ee_ij = diagi(rws_e);  ff_ij = diagi(rws_f);
ef_ij = rws_e + rows*(rws_f-1);

%% Alpha
M.Alpha = Wk./n;

% initialization for Vxx and Y'*Vx
p = P-1;
Xp = cumprod([ones(N,1) x*ones(1,2*p)],2); % x^i => Xp(:,i+1)
Vxx  = zeros(P,P,D);
Vxef = zeros(P,D);
YVx  = zeros(P,D);
gb = zeros(n,2*p+1);  % gamma_mb => gb(:,m+1)
neg = -cumprod(-ones(1,2*p+1));  % neg->(0:2p) => neg(1:2p+1)

% initialization for covariance calculation
zz  = zeros(n,1);
m = 100;  % number of samples
V = zeros(rows);
b   = zeros(n,m,2*p);
e   = zeros(n,m,D,2);
f   = zeros(n,m,D);

% start the M-loop
for k=1:K
  Piik = copy(M.Pik(:,k),lens);

  %% s
  M.S(k) = M.Pik(:,k)'*(M.Eb(:,k).^2 + M.Vb(:,k));
  M.S(k) = M.S(k)./Wk(k);

  %% u
  M.U(k,:) = M.Pik(:,k)'*((M.Ee(:,:,k)-1).^2 + M.Ve(:,:,k));
  M.U(k,:) = M.U(k,:)./Wk(k);

  %% t
  M.T(k,:) = M.Pik(:,k)'*(M.Ef(:,:,k).^2 + M.Vf(:,:,k));
  M.T(k,:) = M.T(k,:)./Wk(k);


  % make power matrices for Eb and Vb
  Ebp = cumprod([ones(n,1) M.Eb(:,k)*ones(1,2*p)],2); % b^i => Ebp(:,i+1)
  Vbp = cumprod([M.Vb(:,k)*ones(1,p)],2);   % Vb^i => Vbp(:,i)
  
  % set up sampling to estimate covariances
  for i=1:n
    Mn = [M.Eb(i,k) M.Ee(i,:,k) M.Ef(i,:,k)];
    V(1)       = M.Vb(i,k);
    V(ee_ij)   = M.Ve(i,:,k);   V(ff_ij)   = M.Vf(i,:,k);
    V(1,rws_e) = M.Vbe(i,:,k);  V(rws_e,1) = M.Vbe(i,:,k)';
    V(1,rws_f) = M.Vbf(i,:,k);  V(rws_f,1) = M.Vbf(i,:,k)';
    V(ef_ij)   = M.Vef(i,:,k);  V(ef_ij-D*(rows-1)) = M.Vef(i,:,k);
    for j=1:MaxPD
     rnd = mvnrnd(Mn,V,m);
      if (isnan(rnd))
        if (j==MaxPD)
          rnd = zeros(m,rows); 
          fprintf('LRM_AA_TT_SH: giving up, posterior covariance is not PD.\n');
          break; 
        end
        V(diagi) = V(diagi)*1.25;   % bias trick, but alleviates probs
      else
        break;
      end
    end
    b(i,:,:) = cumprod(ones(2*p,1)*rnd(:,1)')';  % n-m-2p
    e(i,:,:,1) = rnd(:,rws_e);  % n-m-D-2
    e(i,:,:,2) = e(i,:,:,1).^2;
    f(i,:,:) = rnd(:,rws_f);  % n-m-D
  end
  
  
  % now compute YVx, Vxx, Vxef
  Ee2 = M.Ee(:,:,k).^2;  % we use this below
  for j=0:2*p
    fj = floor(j/2);   % note: floor(1/2)=0, ga(:,2)=gb(:,2)==0
    MjCombo = ones(n,1)*(M.mj(1:fj).*M.Combo(j+1,2*(1:fj)+1));
    NegCombo = ones(N,1)*(neg(1:j+1).*M.Combo(j+1,1:j+1));
    
    % calculate gamma_ma, gamma_mb
    gb(:,j+1) = sum(MjCombo .* Vbp(:,1:fj) .* Ebp(:,j-1:-2:1),2);

    % calculate Gamma_b
    G_b = gb(:,1:j+1);  % n-(j+1)

    for d=1:D
      % D-based covariances
      if (j==0),  Vbpe   = [zz];  
      else
        Vbpe   = [zz M.Vbe(:,d,k) covv2(b(:,:,2:j),e(:,:,d,1))];   % n-(j+1)
      end 
      Vbpe2  = [zz covv2(b(:,:,1:j),e(:,:,d,2))];                  % n-(j+1)
      Vbepf  = [M.Vef(:,d,k) ...
          covv2(b(:,:,1:j).*e(:,:,d(ones(1,j)),1),f(:,:,d))];      % n-(j+1)
      
      % calculate Gammas
      G_be  = G_b.*M.Ee(:,d(ones(1,j+1)),k) + Vbpe;
      G_bef = G_be.*M.Ef(:,d(ones(1,j+1)),k) + Vbepf;
      G_be2 = G_b.*( M.Ve(:,d(ones(1,j+1)),k) + Ee2(:,d(ones(1,j+1))) ) + ...
        Ebp(:,1:j+1).*M.Ve(:,d(ones(1,j+1)),k) + Vbpe2;

      % now, we can calculate the Deltas
      if (j<=p)
        YVx(j+1,d) = Y(:,d)'*(Piik.*sum( ...
          NegCombo.*copy(G_be,lens).*Xp(:,j+1:-1:1) ,2));
        Vxef(j+1,d) = Piik'*sum( ...
          NegCombo.*copy(G_bef,lens).*Xp(:,j+1:-1:1) ,2);
      end
      vxx = Piik'*sum(NegCombo.*copy(G_be2,lens).*Xp(:,j+1:-1:1),2);
      
      % place the results into the matrix form of Vxx
      p1=1; p2=j+1;
      dj = p2-P;  % make sure Vxx is big enough for the whole diagonal
      if (dj>0), p1=p1+dj;  p2=p2-dj; end % if not, then trim it down
      Vxx((P.*([p2:-1:p1]-1) + [p1:p2]) +(d-1)*P^2) = vxx;
    end
  end
  
  % now make Xhat and finish up the M-Step
  Xhat = regmat(x-copy(M.Eb(:,k),lens),P-1);

  for d=1:D
    eXhat = copy(M.Ee(:,d,k),lens)*ones(1,P) .*Xhat;
    PikeXhat = Piik*ones(1,P).*eXhat;
    Ef = copy(M.Ef(:,d,k),lens);
    
    %% Mu
    M.Mu(:,d,k) = (PikeXhat'*eXhat + Vxx(:,:,d)) \ ...
                  (PikeXhat'*(Y(:,d)-Ef) + YVx(:,d) - Vxef(:,d));
    %M.Mu(:,d,k) = (PikeXhat'*eXhat) \ (PikeXhat'*(Y(:,d)-Ef));
    
    % Sigma
    Mu = M.Mu(:,d,k);
    M.Sigma(k,d) = ( Piik'*((Y(:,d)-eXhat*Mu-Ef).^2) ...
      -2*YVx(:,d)'*Mu + Mu'*Vxx(:,:,d)*Mu  ...
      +2*Mu'*Vxef(:,d)+ (M.Pik(:,k).*lens)'*M.Vf(:,d,k) ) ...
      ./ sum(M.Pik(:,k).*lens);
    %M.Sigma(k,d) = (Piik'*((Y(:,d)-eXhat*Mu-Ef).^2)) ...
    %  ./ sum(M.Pik(:,k).*lens);
  end
end
if (M.Options.Sigma.Share)
  M.Sigma = (sum(M.Sigma,2)./D)*ones(1,D);
end
if (M.Options.U.Share)
  M.U = (sum(M.U,2)./D)*ones(1,D);
end
if (M.Options.Tf.Share)
  M.Tf = (sum(M.Tf,2)./D)*ones(1,D);
end

% check for strange circumstances
fs = find(M.S<=0);
fu = find(M.U<=0);  ft = find(M.T<=0);
fsig = find(M.Sigma<=0);
if (~isempty([fs(:); fu(:); ft(:); fsig(:)]))
  fprintf('Some variances are zero\n');
  M.S(fs) = M.Options.minvar;
  M.U(fu) = M.Options.minvar;  M.T(ft) = M.Options.minvar;
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
function M = InitE(M,x,Y,Seq)
[N,D] = size(Y);
n = length(Seq)-1;
K = M.K;
P = M.order+1;

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
M.Ee   = zeros(n,D,K);
M.Ef   = zeros(n,D,K);
M.Vb   = zeros(n,K);
M.Ve   = zeros(n,D,K);
M.Vf   = zeros(n,D,K);
M.Vef  = zeros(n,D,K);
M.Vbe  = zeros(n,D,K);
M.Vbf  = zeros(n,D,K);

M.Eb  = randn(n,K)*.5;
M.Ee  = rand(n,D,K)*2;
M.Ef  = randn(n,D,K)*.5;

M.Vb  = rand(n,K)*2;        % make these positive
M.Ve  = rand(n,D,K)*2;      % make these positive
M.Vf  = rand(n,D,K)*2;      % make these positive
M.Vef = min(M.Ve,M.Vf)./2;  % positive definite

tempVb = permute(M.Vb(:,:,ones(1,D)),[1 3 2]);
M.Vbe = min(tempVb,M.Ve)./2;  % positive definite
M.Vbf = min(tempVb,M.Vf)./2;  % positive definite

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
Ops = SetFieldDef(Ops,'IterLimit',25);
Ops = SetFieldDef(Ops,'NumEMStarts',1);
Ops = SetFieldDef(Ops,'MsgPrefix','');
Ops = SetFieldDef(Ops,'ShowGraphics',0);
Ops = SetFieldDef(Ops,'MsgHnd',[]);
Ops = SetFieldDef(Ops,'MinLen',[]);
Ops = SetFieldDef(Ops,'PropStart',0.8);
Ops = SetFieldDef(Ops,'Sigma',[]);
Ops.Sigma = SetFieldDef(Ops.Sigma,'Share',0);
Ops = SetfieldDef(Ops,'U',[]);
Ops.U = SetFieldDef(Ops.U,'Share',0);
Ops = SetFieldDef(Ops,'Tf',[]);
Ops.Tf = SetFieldDef(Ops.Tf,'Share',0);
Ops = SetFieldDef(Ops,'InitializeWithExamples',1);
%
Ops.SearchOps = optimset('fminsearch');
Ops.SearchOps.Display = 'off';
Ops.SearchOps.MaxFunEvals = 200;
Ops.SearchOps.MaxIter = 200;

function c = covv(x,y)
if (isempty(x) | isempty(y)),  c=[];  return;  end
[m,n,p] = size(x);
c = permute((sum(x.*y,2) - sum(x,2).*sum(y,2)./n)./(n),[1 3 2]);

function c = covv2(x,y)
if (isempty(x) | isempty(y)),  c=[];  return;  end
[m,n,p] = size(x);
y = y(:,:,ones(1,p));
c = squeeze((sum(x.*y,2) - sum(x,2).*sum(y,2)./n)./(n));

function Ops = showmodel(M,trajs)
M.Mu  = permute(M.Mu,[1 3 2]);
M.Ee  = permute(M.Ee,[1 3 2]);
M.Ef  = permute(M.Ef,[1 3 2]);
[trash, M.C] = max(M.Pik,[],2);
Ops = viewmodel(M,trajs,M.Options);

function M = permuteModel(M)
M.Mu  = permute(M.Mu,[1 3 2]);
M.Ee  = permute(M.Ee,[1 3 2]);
M.Ef  = permute(M.Ef,[1 3 2]);
M.Ve  = permute(M.Ve,[1 3 2]);
M.Vf  = permute(M.Vf,[1 3 2]);
M.Vef  = permute(M.Vef,[1 3 2]);
M.Vbe  = permute(M.Vbe,[1 3 2]);
M.Vbf  = permute(M.Vbf,[1 3 2]);

