function [Yhat,Pmod] = lrm_d_ab_pred(M,pt,x,Y,Seq,maxx)
%LRM_D_AB_PRED  Make a partial curve prediction with an LRM_D_AB model.
%
%   This function makes posterior predictions by "predicting" values
%   for all unknown variables. This is in contrast to a likelihood
%   calculation which integrates over (or sums out) all unknown variables.
%   The body of this function is essentially the E-step of the associated
%   cluster model's EM algorithm.
%
%   The main responsibility of this function is to produce partial
%   curve predictions. We take the learned model M and predict the
%   'test' curve point y_hat at x_j using the learned parameters
%   and the partial curve y_i(j-i) (which contains all points up to
%   time j-1). The prediction is calculated in a forward-backward fashion 
%   so that x_j can appear anywhere in the curve.
%
%   As a by-product, this function also returns the posterior model
%   as the second output argument. This model contains all of the 
%   predicted unknown variables (e.g., the membership probabilities)
%   that are required to produce the partial curve prediction.
%   See the code below or the associated EM algorithm for more information.
%
%   [Yhat,PostModel] = LRM_D_AB_PRED(M,pt,X,Y,Seq,['max'])
%    - M       : trained model
%    - pt      : single time point at which to predict y_hat
%    - X,Y,Seq : partial curve in Sequence format (see HELP CCToolbox)
%              : IMPORTANT: length(Seq) MUST equal 2 (i.e., you can only
%              : predict one curve/point with each function call.
%    - max     : see below
%
%   A second calling form is provided that calculates the posterior
%   model for multiple curves simultaneously (i.e., length(Seq)>=2).
%   However, no partial curve prediction is produced in this case and
%   Yhat is returned as empty.
%
%   [[],PostModel] = LRM_D_AB_PRED(M,[],x,Y,Seq,['max'])
%    - M       : trained model
%    - pt      : must equal []
%    - X,Y,Seq : curves in Sequence format (see HELP CCToolbox)
%    - max     : see below
%
%   If you pass the string 'max' as the last argument, then Yhat is
%   calculated from the class w/ maximum membership probability instead
%   of summing across Pik as in the default case.

% Scott Gaffney   10 October 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'lrm_d_ab_pred';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

maxx = cexist('maxx',0);
if (isstr(maxx) & strcmp(maxx,'max'))
  maxx = 1;
else
  maxx = 0;
end

% preprocessing
Mupkd = M.Mu;
M.Mu = permute(M.Mu,[1 3 2]);
[P,D,K] = size(M.Mu);
n = length(Seq)-1;

% Calculate the posterior membership and log-likelihood for the provided
% partial curve information.
Pmod.Ea = ones(n,K);
Pmod.Eb = zeros(n,K);
Pmod.Ef = zeros(n,D,K);
if (isempty(x))
  Pmod.Pik = M.Alpha';   % we are given no curve information so the...
                         % ...posterior membership is just the marginal
  %%%%%%%%%%% Estep
else
  N = Seq(end)-1;
  fun = @postval;
  
  %%%% Calculate posterior mode
  for k=1:K
    r        = M.R(k);
    s        = M.S(k);
    t        = M.T(k,:);
    Mu       = M.Mu(:,:,k);
    sigma    = M.Sigma(k,:);
    SearchOps = M.Options.SearchOps;

    for i=1:n
      indx   = Seq(i):Seq(i+1)-1;
      ni    = length(indx);  
      
      % Find a and b
      pt0 = [1 0]';
      maxpt = fminsearch(fun,pt0,SearchOps,x(indx),Y(indx,:), ...
        Mu,P-1,r,s,t,sigma);
      if (maxpt(1)==0), fprintf('a was zero.\n'); else
        Pmod.Ea(i,k) = maxpt(1);  end;  Pmod.Eb(i,k) = maxpt(2);
      
      % Find f at (a_hat,b_hat)
      Xhat = regmat(M.Ea(i,k)*x(indx)-M.Eb(i,k),P-1);
      Vf  = 1./(ni./sigma + 1./t);
      Pmod.Ef(i,k,:)  = Vf./sigma.*sum(Y(indx,:)-Xhat*Mu);
    end
    
    % Calc Pik
    Pmod.Pik = CalcPik(M,x,Y,Seq);
    s = sum(Pmod.Pik,2);
    Pmod.Lhood_ppt = sum(log(s))./prod(size(Y));
    Pmod.Pik = Pmod.Pik ./ (s*ones(1,K));
    
    % classify sequences
    [trash, Pmod.C] = max(Pmod.Pik,[],2);
  end
end


% Simply return if no prediction is requested
Yhat = [];
if (isempty(pt))
  return;
end

% Generate prediction at pt
if (maxx)
  [trash, k] = max(Pmod.Pik);
  X = regmat(Pmod.Ea(k)*pt-Pmod.Eb(k),P-1);
  Yhat = X*M.Mu(:,:,k) + permute(Pmod.Ef(1,k,:),[1 3 2]);
else
  for d=1:D
    Xk = regmat(Pmod.Ea(1,:)'*pt-Pmod.Eb(1,:)',P-1);
    YhatK = sum(Xk'.*Mupkd(:,:,d)) + Pmod.Ef(1,:,d);
    Yhat(1,d) = Pmod.Pik* YhatK';
  end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% postval 
%%
%%
function val = postval(pt,x,Y,Mu,order,r,s,t,sigma)
a = pt(1);  b = pt(2);
ni = length(x);
D = length(t);
val = 0;
Xhat = regmat(a*x-b,order);
YxMu = Y-Xhat*Mu;
for d=1:D
  iV = eye(ni)/sigma(d) - 1/(ni*sigma(d) + sigma(d)^2/t(d));
  val = val + YxMu(:,d)'*iV*YxMu(:,d);
end
val = val + (a-1)^2/r + b^2/s;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcPik 
%%
function pik = CalcPik(M,x,Y,Seq)
% Numerical integration
NumSamps = 80;
MaxTries = 5;
[N,D] = size(Y);
n = length(Seq)-1;
P = M.order+1;
K = M.K;
Pid = zeros(n,D);
Pik = zeros(n,K);
R = M.R;  S = M.S;  T = M.T;
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
        sigma = M.Sigma(k,d);  t = T(k,d);
        for i=1:n
          indx = Seq(i):Seq(i+1)-1;  ni = length(indx);
          iV = eye(ni)./sigma - 1/(ni*sigma + sigma^2/t);
          Pid(i,d) = mvnormpdf_inv(Y(indx,d)',XMu(indx,d),iV);
        end
      end
      Pik(:,k) = Pik(:,k) + prod(Pid,2);
    end
  end
  pik = (Pik./TotalSamps) .* (ones(n,1)*M.Alpha'); % we keep Pik for next try!
  if (all(sum(pik,2))), break; end

  % we have detected some zeros, try again?
  if (tries==MaxTries)
    fprintf('Integration failed, using realmin*1e100 instead.\n');
    zero = find(sum(pik,2)==0);
    pik(zero,:) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
    break;
  else
    fprintf('Zero membership detected, trying integration again: %d\n',tries);
    tries = tries+1;
    S = 1.5*S;  % biased, but gets over some tricky integrations
    R = 1.5*R;
  end
end
