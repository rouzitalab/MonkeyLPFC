function [Yhat,Pmod] = lrm_pred(M,pt,x,Y,Seq,maxx)
%LRM_PRED  Make a partial curve prediction with an LRM model.
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
%   [Yhat,PostModel] = LRM_PRED(M,pt,X,Y,Seq,['max'])
%    - M       : trained model
%    - pt      : single time point at which to predict y_hat
%    - X,Y,Seq : partial curve in Sequence format; See also CCTOOLBOX
%              : IMPORTANT: length(Seq) must equal 2, i.e., you can only
%              : predict one curve/point with each function call.
%    - max     : see below
%
%   A second calling form is provided that calculates the posterior
%   model for multiple curves simultaneously (i.e., length(Seq)>=2).
%   However, no partial curve prediction is produced in this case and
%   Yhat is returned as empty.
%
%   [[],PostModel] = LRM_PRED(M,[],x,Y,Seq,['max'])
%    - M       : trained model
%    - pt      : must equal []
%    - X,Y,Seq : curves in Sequence format; See also CCTOOLBOX
%    - max     : see below
%
%   If you pass the string 'max' as the last argument, then Yhat is
%   calculated from the class w/ maximum membership probability instead
%   of summing across Pik as in the default case.

% Scott Gaffney   10 October 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'lrm_pred';
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
Mu = permute(M.Mu,[1 3 2]);
[P,D,K] = size(Mu);
likefunc = listmodels(mfilename,'like');

% Calculate the posterior membership and log-likelihood for the provided
% partial curve information.
if (isempty(x))
  Pmod.Pik = M.Alpha';  % we are given no curve information so the...
else                    % ...posterior membership is just the marginal
  [Pmod.Lhood_ppt, other] = feval(likefunc,M,x,Y,Seq);
  Pmod.Pik = other.Pik;
  Pmod.C = other.C;
end

% Simply return if no prediction is requested
Yhat = [];
if (isempty(pt))
  return;
end

% Generate prediction at pt (pt is scalar)
X = regmat(pt,P-1);
if (maxx)
  [trash, k] = max(Pmod.Pik);
  Yhat = X*Mu(:,:,k);
else
  for d=1:D
    YhatK = X*Mupkd(:,:,d);
    Yhat(1,d) = Pmod.Pik* YhatK';
  end
end
