function [m_hat,sigma,B] = kernreg_mv(X,Y,Q,kernel,h,W)
%KERNREG_MV  Perform multivariate (output) kernel regression.
%   KERNREG_MV(X,Y,Q,kernel,H,W) is a vector of predicted Y-values at the
%   query points contained in the vector Q by way of kernel regression.
%   X and Y are matrices containing the (x,y) data pairs that will be
%   used in performing weighted least squares.  Y is N by D and X is N by P,
%   where N is the number of data points, D is the dimension of the output
%   vector, and P-1 is the order of the model (this assumes that X(:,1) is
%   a vector of ones. W is an optional weight vector that can be specified; 
%   if W~=[], its elements will be multiplied by the corresponding kernel 
%   weights for each data point. H optionally specifies the bandwidth used 
%   with the kernels. The type of kernel used in the regression can be 
%   specified by the string argument kernel as follows.
%
%     'normal'      : gaussian pdf <default>
%     'uniform'     : uniform pdf
%     'pass'        : pass through (the kernel weights are in the 
%                     matrix W already)
%
%  [PRED,SIGMA,B] = KERNREG_MV(...) behaves as above and in addition returns
%  the estimated covariance matrices for the error terms in SIGMA
%  and the regression coefficients as column vectors in the 
%  matrix B. SIGMA is D by D by N, and B is P by D.

% Scott J. Gaffney   20 January 1999
% DataLab @UCI
% Department of Information and Computer Science
% University of California, Irvine, USA.

PROGNAME = 'kernreg_mv';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

NORMAL_BANDWIDTH     = .3;
UNIFORM_BANDWIDTH    = .3;

%% Preprocessing
if (nargin < 3)
  error('kernreg_mv: not enough arguments were specified to the function.');
end
Q=Q(:);
sizeX = size(X);
sizeY = size(Y);
N = sizeX(1);
numQry = length(Q);
Order = sizeX(2)-1;
h = cexist('h',[]);
kernel = cexist('kernel','normal');
W = cexist('W',ones(N,1));

if (sizeX(1) ~= sizeY(1))
  error('kernreg_mv: size(X,1) must equal size(Y,1).\n');
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
      K(:,i) = normpdf((X(:,2)-Q(i))./h);
    end
end


%**********************************************************************
% Perform Kernel Regression
%
pointVector = ones(1,Order);
if (prod(size(W)) == length(W))   % W is a vector?
  for i=1:numQry
    [B(:,:,i),sigma(:,:,i)] = wls(X,Y,W.*K(:,i));
    m_hat(i,:) = cumprod([1 pointVector*Q(i)]) * B(:,:,i);
  end
else
  for i=1:numQry
    [B(:,:,i),sigma(:,:,i)] = wls(X,Y,W(:,i).*K(:,i));
    m_hat(i,:) = cumprod([1 pointVector*Q(i)]) * B(:,:,i);
  end
end

