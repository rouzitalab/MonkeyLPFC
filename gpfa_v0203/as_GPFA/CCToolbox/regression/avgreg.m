function [m_hat,sigma] = avgreg(X,Y,Q,W)
%AVGREG  Perform regression by taking pointwise averages.
%   AVGREG(X,Y,Q,W) is a vector of predicted Y-values at the
%   query points contained in the vector Q by way of pointwise averaging.
%   That is, the predicted value at Q(1) is mean(Y(find(X==Q(1))));
%
%   X and Y are vectors containing the (x,y) data pairs. 
%   W is an optional weight vector that can be specified.
%   if W~=[], its elements will be used to compute pointwise
%   weighted averages.
%
%   [PRED,SIGMA] = AVGREG(...) behaves as above and in addition returns
%   the estimated standard deviation at each query point about the 
%   predicted Y-value.

% Scott J. Gaffney   27 January 1999
% DataLab @UCI
% Department of Information and Computer Science
% University of California, Irvine, USA.

%% Preprocessing
if (nargin < 3)
  error('avgreg: not enough arguments were specified to the function.');
end
Y=Y(:);
Q=Q(:);
N = size(X,1);
numQry = length(Q);

if (N ~= length(Y))
  error('avgreg: size(X,1) must equal size(Y,1).');
end

if (exist('W','var') ~= 1 | isempty(W))
  W = ones(N,1);
else
  W=W(:);
  if (length(W) ~= N)
    error('avgreg: W must be a vector of length size(X,1).');
  end
end


wY = Y .* W;
for i=1:numQry
  tmp = wY(find(X==Q(i)));
  m_hat(i) = mean(tmp);
  sigma(i) = var(tmp);
end

sigma = sigma(:);
m_hat = m_hat(:);
