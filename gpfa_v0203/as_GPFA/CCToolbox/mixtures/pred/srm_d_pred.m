function [Yhat,Pmod] = srm_d_pred(M,pt,x,Y,Seq,maxx)
%SRM_D_PRED  Make a partial curve prediction with an SRM_D model.
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
%   [Yhat,PostModel] = SRM_D_PRED(M,pt,X,Y,Seq,['max'])
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
%   [[],PostModel] = SRM_D_PRED(M,[],x,Y,Seq,['max'])
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

PROGNAME = 'srm_d_pred';
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
Pmod.Ef = zeros(n,K,D);
if (isempty(x))
  Pmod.Pik = M.Alpha';   % we are given no curve information so the...
                         % ...posterior membership is just the marginal
else
  X = bsplinebasis(M.knots,M.order,x);
  Pikd = zeros(n,K,D);
  for k=1:K
    for d=1:D
      Mu       = M.Mu(:,d,k);
      sigma    = M.Sigma(k,d);
      s        = M.S(k,d);
      for i=1:n
        indx   = Seq(i):Seq(i+1)-1;
        ni     = length(indx);
        XMu    = X(indx,:)*Mu;
        iS     = eye(ni)/sigma - 1/(ni*sigma + sigma^2/s);
        
        Pikd(i,k,d) = mvnormpdf_inv(Y(indx,d)',XMu',iS);
        Pmod.Ef(i,k,d)  = s/(ni*s+sigma)*sum(Y(indx,d)-XMu);
      end
    end
  end
  Pmod.Pik = prod(Pikd,3).* (ones(n,1)*M.Alpha');
  s = sum(Pmod.Pik,2);
  if (~all(s))
    fprintf([PROGNAME, ': log(0) detected, using log(K*realmin*1e100).\n']);
    zero = find(s==0);
    Pmod.Pik(zero,:) = realmin*1e100*(ones(length(zero),1)*M.Alpha');
    s(zero) = sum(Pmod.Pik(zero,:),2);
  end
  Pmod.Lhood_ppt = sum(log(s))./prod(size(Y));
  Pmod.Pik = Pmod.Pik ./ (s*ones(1,K));
  [trash, Pmod.C] = max(Pmod.Pik,[],2);
end

% Simply return if no prediction is requested
Yhat = [];
if (isempty(pt))
  return;
end

% Generate prediction at pt
X = bsplinebasis(M.knots,M.order,pt);
if (maxx)
  [trash, k] = max(Pmod.Pik);
  Yhat = X*M.Mu(:,:,k) + permute(Pmod.Ef(1,k,:),[1 3 2]);
else
  for d=1:D
    YhatK = X*Mupkd(:,:,d) + Pmod.Ef(1,:,d);
    Yhat(1,d) = Pmod.Pik* YhatK';
  end
end





