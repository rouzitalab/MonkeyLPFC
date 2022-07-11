function [Yhat,Pmod] = srm_ab_pred(M,pt,x,Y,Seq,maxx)
%SRM_AB_PRED  Make a partial curve prediction with an SRM_AB model.
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
%   [Yhat,PostModel] = SRM_AB_PRED(M,pt,X,Y,Seq,['max'])
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
%   [[],PostModel] = SRM_AB_PRED(M,[],x,Y,Seq,['max'])
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

PROGNAME = 'srm_ab_pred';
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
if (isempty(x))
  Pmod.Pik = M.Alpha';   % we are given no curve information so the...
                         % ...posterior membership is just the marginal
%%%%%%%%%%% Estep
else
  N = Seq(end)-1;
  fun = @postval;
  SearchOps = M.Options.SearchOps;
  
  %%%% Calculate posterior mode
  for k=1:K
    for i=1:n
      indx = Seq(i):Seq(i+1)-1;
      pt0 = [1;0];
      %pt0 = [(1+M.Ea(i,k))/2 M.Eb(i,k)/2]';
      maxpt = fminsearch(fun,pt0,SearchOps,x(indx),Y(indx,:),M.Mu(:,:,k), ...
        M.knots,M.order,M.R(k),M.S(k),M.Sigma(k,:));
      Pmod.Ea(i,k) = maxpt(1);
      Pmod.Eb(i,k) = maxpt(2);
    end
  end

  % Calc Pik
  [Pmod.Pik,scale] = CalcPik(M,x,Y,Seq);
  s = sum(Pmod.Pik,2);
  Pmod.Lhood_ppt = (sum(log(s)) + N*log(scale))./prod(size(Y));
  Pmod.Pik = Pmod.Pik ./ (s*ones(1,K));

  % classify sequences
  [trash, Pmod.C] = max(Pmod.Pik,[],2);
end


% Simply return if no prediction is requested
Yhat = [];
if (isempty(pt))
  return;
end

% Generate prediction at pt
if (maxx)
  [trash, k] = max(Pmod.Pik);
  X = bsplinebasis(M.knots,M.order,Pmod.Ea(k)*pt-Pmod.Eb(k));
  Yhat = X*M.Mu(:,:,k);
else
  for d=1:D
    Xk = bsplinebasis(M.knots,M.order,Pmod.Ea(1,:)'*pt-Pmod.Eb(1,:)');
    YhatK = sum(Xk'.*Mupkd(:,:,d));
    Yhat(1,d) = Pmod.Pik* YhatK';
  end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% postval 
%%
%%
function val = postval(pt,x,Y,Mu,knots,order,r,s,sigma)
a = pt(1);  b = pt(2);
Xhat = bsplinebasis(knots,order,a*x-b);
val = sum(sum((Y-Xhat*Mu).^2)./sigma) + (a-1)^2/r + b^2/s;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CalcPik 
%%
function [Pik,scale] = CalcPik(M,x,Y,Seq)
% Numerical integration
NumSamps = 80;
MaxTries = 5;
[N,D] = size(Y);
n = length(Seq)-1;
K = M.K;
mlen = max(diff(Seq));
Piid = zeros(N,D);
Piimk = zeros(N,NumSamps,K);
Pik = zeros(n,K);
S = M.S;  R = M.R;
TotalSamps = 0;
tries = 1;

while (1)
  TotalSamps = TotalSamps + NumSamps;
  Pik(:) = 0;

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
  scale = mean(mean(mean(Piimk)));
  Piimk_scl = Piimk./scale;
  for k=1:K
    for j=1:TotalSamps
      Pik(:,k) = Pik(:,k) + sprod(Piimk_scl(:,j,k),Seq,mlen);
    end
  end
  clear Piimk_scl;
  Pik = (Pik./TotalSamps) .* (ones(n,1)*M.Alpha');
  if (all(sum(Pik,2))), break; end

  % we have detected some zeros, try again?
  if (tries==MaxTries)
    fprintf(['srm_tt_sh_pred: Integration failed, using realmin*1e100 ',...
      'instead.\n']);
    zero = find(sum(Pik,2)==0);
    Pik(zero,:) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
    break;
  else
    fprintf(['srm_tt_sh_pred: Zero membership detected, trying ', ...
        'integration again: %d\n'],tries);
    tries = tries+1;
    S = 1.25*S;  % biased, but gets over some tricky integrations
    R = 1.25*R;  % biased, but gets over some tricky integrations
    Piimk = [zeros(N,NumSamps,K) Piimk]; % save current values
  end
end
