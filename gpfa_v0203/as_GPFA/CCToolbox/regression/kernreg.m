function [m_hat,sigma,B] = kernreg(X,Y,Q,kernel,h,W,NoNorm)
%KERNREG  Perform kernel regression.
%   KERNREG(X,Y,Q,kernel,H,W) is a vector of predicted Y-values at the
%   query points contained in the vector Q by way of kernel regression.
%   X and Y are vectors containing the (x,y) data pairs that will be
%   used in performing weighted least squares.  W is an optional weight
%   matrix that can be specified; if W~=[], its elements will be
%   multiplied by the corresponding kernel weights for each data point.
%   H optionally specifies the bandwidth used with the kernels.
%   The type of kernel used in the regression can be specified by
%   the string argument kernel as follows.
%
%     'normal'      : gaussian pdf <default>
%     'uniform'     : uniform pdf
%     'pass'        : pass through (the kernel weights are in the 
%                     matrix W already)
%
%  [PRED,SIGMA,B] = KERNREG(...,NoNorm) behaves as above and in addition returns
%  the estimated variance about the regression function in 
%  sigma, and the regression coefficients as column vectors in the 
%  matrix B. If NoNorm does not exist or if it is not equal to 1, SIGMA
%  is normalized as in WLS.

% Scott J. Gaffney   20 January 1999
% DataLab @UCI
% Department of Information and Computer Science
% University of California, Irvine, USA.
%
% (SJG) 18 December 2000
% Added the ability to have an N-by-N W matrix.
%
% (SJG) 20 December 2000
% Added the ability to pass through the kernel weights in W.
%

NORMAL_BANDWIDTH     = .3;
UNIFORM_BANDWIDTH    = .3;

%% Preprocessing
if (nargin < 3)
  error('kernreg: not enough arguments were specified to the function.\n');
end
Q=Q(:);
Y=Y(:);
sizeX = size(X);
sizeY = size(Y);
N = sizeX(1);
numQry = length(Q);
Order = sizeX(2)-1;  % there should be a column of ones, thus subtract 1
h = cexist('h',[]);
kernel = cexist('kernel','normal');
W = cexist('W',ones(N,1));

if (sizeX(1) ~= sizeY(1))
  error('kernreg: size(X,1) must equal size(Y,1).\n');
end

if (exist('NoNorm','var') ~= 1 | NoNorm ~= 1)
  NoNorm = 0;
else
  NoNorm = 1;
end


%% Calculate the Kernel Weights
switch(kernel)
  case 'uniform',
    if (isempty(h)), h = UNIFORM_BANDWIDTH; end;
    for i=1:numQry
      K(:,i) = unipdf(abs(X(:,2)-Q(i)),0,h);
    end
  case 'pass'
    K = ones(N,numQry);
  case 'normal',
    if (isempty(h)), h = NORMAL_BANDWIDTH; end;
    for i=1:numQry
      K(:,i) = normpdf((X(:,2)-Q(i))./h);  % ./h, no need for this
    end
end




%**********************************************************************
% Perform Kernel Regression
%
pointVector = ones(1,Order);
if (prod(size(W)) == length(W))   % W is a vector?
  for i=1:numQry
    [B(:,i),sigma(i)] = wls(X,Y,W.*K(:,i),NoNorm);
    m_hat(i) = cumprod([1 pointVector*Q(i)]) * B(:,i);
  end
else
  for i=1:numQry
    [B(:,i),sigma(i)] = wls(X,Y,W(:,i).*K(:,i),NoNorm);
    m_hat(i) = cumprod([1 pointVector*Q(i)]) * B(:,i);
  end
end

sigma = sigma(:);
m_hat = m_hat(:);
