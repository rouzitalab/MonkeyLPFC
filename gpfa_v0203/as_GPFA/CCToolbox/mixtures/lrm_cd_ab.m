function M = lrm_cd_ab(trajs,K,order,Ops)
%LRM_CD_AB  LRM with transformation y=c[ax+b]+d
%
%   Model = LRM_CD_AB(Trajs,K,order,[Options])
%    - Trajs : 'Trajs' structure (see HELP CCToolbox)
%    - K     : number of clusters to find
%    - order : order of polynomial regression
%
%   DefOps = LRM_CD_AB('options) returns the default options.

% Scott Gaffney   1 Aug 2003
% Department of Information and Computer Science
% University of California, Irvine
%
% Changes
% -----------------

PROGNAME = 'lrm_cd_ab';
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
M.NumIndParams = (K-1) + K*P*D + 3*K*D + 2*K;% alpha, mu, sigma,u,t, r,s
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
function val = b_y(pt,x,Y,Mu,order,r,s,u,t,sigma)
a = pt(1);  b = pt(2);
ni = length(x);
D = length(t);
val = 0;
Xhat = regmat(a*x-b,order);
XMu = Xhat*Mu;
YxMu = Y-XMu;
I = eye(ni);
for d=1:D
  iR = I/sigma(d) - XMu(:,d)*XMu(:,d)'/(sigma(d)^2/u(d) ...
    + sigma(d)*(XMu(:,d)'*XMu(:,d)));
  siR = sum(iR);
  iV = iR - sum(iR,2)*siR/(1/t(d) + sum(siR));
  val = val + YxMu(:,d)'*iV*YxMu(:,d);
end
val = val + (a-1)^2/r + b^2/s;


function val = post(pt,x,Y,Mu,order,r,s,u,t,sigma)
D = length(t);
a = pt(1);  b = pt(2);  e = pt(3:3+D-1);  f = pt(3+D:end);
ni = length(x);
val = 0;
XMu = regmat(a*x-b,order)*Mu;
ee = (e-1).^2./u;  ff = f.^2./t;
for d=1:D
  val = val + sum((Y(:,d)-e(d)*XMu(:,d)-f(d)).^2)./sigma(d) + ee(d) + ff(d);
end
val = val + (a-1)^2/r + b^2/s;



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
rows  = 2+2*D;         diagi = [1:rows+1:rows^2];
rws_e = 3:3+D-1;       rws_f = 3+D:3+D+D-1;
aa_ij = diagi(1);      bb_ij = diagi(2);
ee_ij = diagi(rws_e);  ff_ij = diagi(rws_f);
ab_ij = 2;             ef_ij = rws_e + rows*(rws_f-1);
ae_ij = rws_e;         af_ij = rws_f;
be_ij = rws_e+rows;    bf_ij = rws_f+rows;


for k=1:K
  r        = M.R(k);
  s        = M.S(k);
  u        = M.U(k,:);
  t        = M.T(k,:);
  Mu       = M.Mu(:,:,k);
  sigma    = M.Sigma(k,:);
  for i=1:n
    indx   = Seq(i):Seq(i+1)-1;
    ni    = length(indx);  
    
    a0 = M.Ea(i,k) - (1-M.Options.PropStart)*(M.Ea(i,k)-1);
    b0 = M.Eb(i,k) * M.Options.PropStart;
    pt0 = [a0 b0]';
    maxpt = fminsearch(fun,pt0,SearchOps,x(indx),Y(indx,:), ...
      Mu,P-1,r,s,u,t,sigma);
    if (maxpt(1)==0), fprintf('a was zero.\n'); else
    M.Ea(i,k) = maxpt(1);  end;  M.Eb(i,k) = maxpt(2);

    % Find e and f at (a_hat,b_hat)
    Xhat = regmat(M.Ea(i,k)*x(indx)-M.Eb(i,k),P-1);
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
    
    % Now we will calc the inverse Information at Ea,Eb and Ee,Ef.
    Dx = ones(ni,1)*(0:P-1);
    Dx(:,3:end) = Dx(:,3:end) .* Xhat(:,2:end-1);
    D2x = ones(ni,1)*[0 0 2:P-1];
    D2x(:,3:end) = D2x(:,3:end) .* Dx(:,2:end-1);
    
    % ...then make repeat calculations
    Mue = Mu.*(ones(P,1)*M.Ee(i,:,k));
    f = ones(ni,1)*M.Ef(i,:,k);
    Xmu     = Xhat*Mu; % ni x D
    Dxmu    = Dx*Mu;  % ni x D
    Dxmux   = Dxmu.*x(indx,ones(1,D));  % ni x D
    DxMu    = Dx*Mue;  % ni x D
    DxMux   = DxMu.*x(indx,ones(1,D));  % ni x D
    D2xMu   = D2x*Mue;  % ni x D
    D2xMux  = D2xMu.*x(indx,ones(1,D));  % ni x D
    D2xMux2 = D2xMux.*x(indx,ones(1,D));  % ni x D
    YxMu = Y(indx,:) - Xhat*Mue - f;    % ni x D

    % ...and now we can calculate the information.
    Ia  = sum(sum(YxMu.*D2xMux2)./sigma) - sum(sum(DxMux.^2)./sigma) - 1./r;
    Ib  = sum(sum(YxMu.*D2xMu)./sigma)   - sum(sum(DxMu.^2)./sigma)  - 1./s;
    Iab = -sum(sum(YxMu.*D2xMux)./sigma) + sum(sum(DxMu.*DxMux)./sigma);
    Iae = (sum(YxMu.*Dxmux) - sum(Xmu.*DxMux))./sigma;
    Ibe = -(sum(YxMu.*Dxmu) - sum(Xmu.*Dxmu))./sigma;
    Iaf = -sum(DxMux)./sigma;
    Ibf = sum(DxMu)./sigma;
    Ie  = -sum(Xmu.^2)./sigma - 1./u;
    If  = -ni./sigma - 1./t;
    Ief = -sum(Xmu)./sigma;
    
    % the variance is the inverse of the negative information
    V = diag([Ia Ib Ie If]);
    V(1,2)     = Iab;  V(2,1)     = Iab;
    V(1,rws_e) = Iae;  V(rws_e,1) = Iae';
    V(2,rws_e) = Ibe;  V(rws_e,2) = Ibe';
    V(1,rws_f) = Iaf;  V(rws_f,1) = Iaf';
    V(2,rws_f) = Ibf;  V(rws_f,2) = Ibf';
    V(ef_ij)   = Ief;  V(ef_ij-D*(rows-1)) = Ief;

    V = inv(-V);
    M.Va(i,k)    = V(aa_ij);
    M.Vb(i,k)    = V(bb_ij);
    M.Ve(i,:,k)  = V(ee_ij);
    M.Vf(i,:,k)  = V(ff_ij);
    M.Vab(i,k)   = V(ab_ij);
    M.Vef(i,:,k) = V(ef_ij);
    M.Vae(i,:,k) = V(ae_ij);
    M.Vbe(i,:,k) = V(be_ij);
    M.Vaf(i,:,k) = V(af_ij);
    M.Vbf(i,:,k) = V(bf_ij);
  end
end
fa = find(M.Va<=0);    fb = find(M.Vb<=0);
fe = find(M.Ve<=0);    ff = find(M.Vf<=0);
if (~isempty([fa(:); fb(:); fe(:); ff(:)]))
  fprintf('A negative variance was detected....\n');
  M.Va(fa) = -M.Va(fa);  M.Vb(fb) = -M.Vb(fb);
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
R = M.R;  S = M.S;  U = M.U;  T = M.T;
TotalSamps = 0;
tries = 1;

while (1)
  TotalSamps = TotalSamps + NumSamps;

  % calculate the density at sampled points
  for k=1:K
    a = randn(NumSamps,1).*sqrt(R(k)) + 1;  % sample from N(1,r)
    b = randn(NumSamps,1).*sqrt(S(k));   % sample from N(0,s)
    for j=1:NumSamps
      Xhat = regmat(a(j)*x-b(j),P-1);
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
    S = 1.5*S;  % biased, but gets over some tricky integrations
    R = 1.5*T;
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
    fprintf('The log-likelihood is equal to NaN.\n');
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
rows  = 2+2*D;         diagi = [1:rows+1:rows^2];
rws_e = 3:3+D-1;       rws_f = 3+D:3+D+D-1;
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
ga = zeros(n,2*p+1);  % gamma_ma => ga(:,m+1)
gb = zeros(n,2*p+1);  % gamma_mb => gb(:,m+1)
neg = -cumprod(-ones(1,2*p+1));  % neg->(0:2p) => neg(1:2p+1)

% initialization for covariance calculation
zz  = zeros(n,1);
m = 100;  % number of samples
V = zeros(rows);
a   = zeros(n,m,2*p);
b   = zeros(n,m,2*p);
e   = zeros(n,m,D,2);
f   = zeros(n,m,D);

% start the M-loop
for k=1:K
  Piik = copy(M.Pik(:,k),lens);

  %% r
  M.R(k) = M.Pik(:,k)'*((M.Ea(:,k)-1).^2 + M.Va(:,k));
  M.R(k) = M.R(k)./Wk(k);

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
  Eap = cumprod([ones(n,1) M.Ea(:,k)*ones(1,2*p)],2); % a^i => Eap(:,i+1)
  Ebp = cumprod([ones(n,1) M.Eb(:,k)*ones(1,2*p)],2); % b^i => Ebp(:,i+1)
  Vap  = cumprod([M.Va(:,k)*ones(1,p)],2);   % Va^i => Vap(:,i)
  Vbp  = cumprod([M.Vb(:,k)*ones(1,p)],2);   % Vb^i => Vbp(:,i)
  
  % set up sampling to estimate cov(a^(p-i),b^i), cov(a^(p-i)*b^i,e^{1|2}), and
  % cov(a^(p-i)*b^i*e,f)
if (1)  
  for i=1:n
    Mn = [M.Ea(i,k) M.Eb(i,k) M.Ee(i,:,k) M.Ef(i,:,k)];
    V(1)       = M.Va(i,k);     V(2,2)     = M.Vb(i,k);
    V(ee_ij)   = M.Ve(i,:,k);   V(ff_ij)   = M.Vf(i,:,k);
    V(1,2)     = M.Vab(i,k);    V(2,1)     = M.Vab(i,k);
    V(1,rws_e) = M.Vae(i,:,k);  V(rws_e,1) = M.Vae(i,:,k)';
    V(2,rws_e) = M.Vbe(i,:,k);  V(rws_e,2) = M.Vbe(i,:,k)';
    V(1,rws_f) = M.Vaf(i,:,k);  V(rws_f,1) = M.Vaf(i,:,k)';
    V(2,rws_f) = M.Vbf(i,:,k);  V(rws_f,2) = M.Vbf(i,:,k)';
    V(ef_ij)   = M.Vef(i,:,k);  V(ef_ij-D*(rows-1)) = M.Vef(i,:,k);
    for j=1:MaxPD
      rnd = mvnrnd(Mn,V,m);
      if (isnan(rnd))
        if (j==MaxPD)
          rnd = zeros(m,rows); 
          %fprintf('LRM_AA_TA_SH: giving up, posterior covariance is not PD.\n');
          break; 
        end
        V(diagi) = V(diagi)*1.25;   % bias trick, but alleviates probs
      else
        break;
      end
    end
    a(i,:,:) = cumprod(ones(2*p,1)*rnd(:,1)')';  % n-m-2p
    b(i,:,:) = cumprod(ones(2*p,1)*rnd(:,2)')';  % n-m-2p
    e(i,:,:,1) = rnd(:,rws_e);  % n-m-D-2
    e(i,:,:,2) = e(i,:,:,1).^2;
    f(i,:,:) = rnd(:,rws_f);  % n-m-D
  end
  
  
  % now compute YVx, Vxx, Vxef
  Ee2 = M.Ee(:,:,k).^2;
  for j=0:2*p
    fj = floor(j/2);   % note: floor(1/2)=0, ga(:,2)=gb(:,2)==0
    MjCombo = ones(n,1)*(M.mj(1:fj).*M.Combo(j+1,2*(1:fj)+1));
    NegCombo = ones(N,1)*(neg(1:j+1).*M.Combo(j+1,1:j+1));
    
    % calculate gamma_ma, gamma_mb
    ga(:,j+1) = sum(MjCombo .* Vap(:,1:fj) .* Eap(:,j-1:-2:1),2);
    gb(:,j+1) = sum(MjCombo .* Vbp(:,1:fj) .* Ebp(:,j-1:-2:1),2);

    % calculate cov(a^(j-r),b^r)
    if     (j==0), Vapbp = [zz];
    elseif (j==2), Vapbp = [zz M.Vab(:,k) zz];  % n-(j+1)
    else           Vapbp = [zz covv(a(:,:,j-1:-1:1),b(:,:,1:j-1)) zz];
    end

    % calculate Gamma_ab
    G_ab = Eap(:,j+1:-1:1).*gb(:,1:j+1) + Ebp(:,1:j+1).*ga(:,j+1:-1:1) + ...
      ga(:,j+1:-1:1).*gb(:,1:j+1) + Vapbp;  % n-(j+1)

    if (j>0)
      ab = a(:,:,j);  ab(:,:,j+1) = b(:,:,j);
      ab(:,:,2:j) = a(:,:,j-1:-1:1).*b(:,:,1:j-1);  % n-m-(j+1)
    end
    for d=1:D
      % cov(a^(j-r)*b^r,e), cov(a^(j-r)*b^r,e^2), cov(a^(j-r)*b^r*e,f)
      if (j==0)
        Vabpe = zz;  Vabpe2 = zz;  Vabepf = M.Vef(:,d,k);
      else
        if (j==1)
          Vabpe = [M.Vae(:,d,k) M.Vbe(:,d,k)];
        else
          Vabpe = covv2(ab,e(:,:,d,1));
        end
        Vabpe2  = covv2(ab,e(:,:,d,2));                         % n-(j+1)
        Vabepf  = covv2(ab.*e(:,:,d(ones(1,j+1)),1),f(:,:,d));  % n-(j+1)
      end
      
      % calculate Gamma_abe, Gamma_abe^2, and Gamma_abef
      G_abe  = G_ab.*M.Ee(:,d(ones(1,j+1)),k) + Vabpe;
      G_abef = G_abe.*M.Ef(:,d(ones(1,j+1)),k) + Vabepf;
      G_abe2 = G_ab.*( M.Ve(:,d(ones(1,j+1)),k) + Ee2(:,d(ones(1,j+1))) ) + ...
        Eap(:,j+1:-1:1).*Ebp(:,1:j+1).*M.Ve(:,d(ones(1,j+1)),k) + Vabpe2;

      % now, we can calculate the Deltas
      if (j<=p)
        YVx(j+1,d) = Y(:,d)'*(Piik.*sum( ...
          NegCombo.*copy(G_abe,lens).*Xp(:,j+1:-1:1) ,2));
        Vxef(j+1,d) = Piik'*sum( ...
          NegCombo.*copy(G_abef,lens).*Xp(:,j+1:-1:1) ,2);
      end
      vxx = Piik'*sum( NegCombo.*copy(G_abe2,lens).*Xp(:,j+1:-1:1) ,2);
      
      % place the results into the matrix form of Vxx
      p1=1; p2=j+1;
      dj = p2-P;  % make sure Vxx is big enough for the whole diagonal
      if (dj>0), p1=p1+dj;  p2=p2-dj; end % if not, then trim it down
      Vxx((P.*([p2:-1:p1]-1) + [p1:p2])+(d-1)*P^2) = vxx;
    end
  end
end  % if (1), a debugging thing

  % now make Xhat and finish up the M-Step
  Xhat = regmat(copy(M.Ea(:,k),lens).*x-copy(M.Eb(:,k),lens),P-1);

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
    -2*YVx(:,d)'*Mu+ Mu'*Vxx(:,:,d)*Mu  ...
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
fr = find(M.R<=0);  fs = find(M.S<=0);
fu = find(M.U<=0);  ft = find(M.T<=0);
fsig = find(M.Sigma<=0);
if (~isempty([fr(:); fs(:); fu(:); ft(:); fsig(:)]))
  fprintf('Some variances are zero\n');
  M.R(fr) = M.Options.minvar;  M.S(fs) = M.Options.minvar;
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
M.Ea   = zeros(n,K);
M.Eb   = zeros(n,K);
M.Ee   = zeros(n,D,K);
M.Ef   = zeros(n,D,K);
M.Va   = zeros(n,K);
M.Vb   = zeros(n,K);
M.Ve   = zeros(n,D,K);
M.Vf   = zeros(n,D,K);
M.Vab  = zeros(n,K);
M.Vef  = zeros(n,D,K);
M.Vae  = zeros(n,D,K);
M.Vaf  = zeros(n,D,K);
M.Vbe  = zeros(n,D,K);
M.Vbf  = zeros(n,D,K);

M.Ea  = rand(n,K)*2;
M.Eb  = randn(n,K)*.5;
M.Ee  = rand(n,D,K)*2;
M.Ef  = randn(n,D,K)*.5;

M.Va  = rand(n,K)*2;        % make these positive
M.Vb  = rand(n,K)*2;        % make these positive
M.Ve  = rand(n,D,K)*2;      % make these positive
M.Vf  = rand(n,D,K)*2;      % make these positive
M.Vab = min(M.Va,M.Vb)./2;  % positive definite
M.Vef = min(M.Ve,M.Vf)./2;  % positive definite

tempVa = permute(M.Va(:,:,ones(1,D)),[1 3 2]);
tempVb = permute(M.Vb(:,:,ones(1,D)),[1 3 2]);
M.Vae = min(tempVa,M.Ve)./2;  % positive definite
M.Vaf = min(tempVa,M.Vf)./2;  % positive definite
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
Ops = SetFieldDef(Ops,'Sigma',[]);
Ops.Sigma = SetFieldDef(Ops.Sigma,'Share',0);
Ops = SetfieldDef(Ops,'Tf',[]);
Ops.Tf = SetFieldDef(Ops.Tf,'Share',0);
Ops = SetFieldDef(Ops,'U',[]);
Ops.U = SetFieldDef(Ops.U,'Share',0);
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

function c = covv2(x,y)
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
M.Vae  = permute(M.Vae,[1 3 2]);
M.Vaf  = permute(M.Vaf,[1 3 2]);
M.Vbe  = permute(M.Vbe,[1 3 2]);
M.Vbf  = permute(M.Vbf,[1 3 2]);

