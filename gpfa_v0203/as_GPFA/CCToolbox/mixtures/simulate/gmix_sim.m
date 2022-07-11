function [X,M,MLMU,MLSIGMA] = gmix_sim(M,N)
%GMIX_SIM  Generates multivariate normal data from multiple classes.
%   [X,M,MLMU,MLSIGMA] = GMIX_SIM(M,N) generates a data set 
%   with N rows of data points, generated from K gaussian 
%   distributions according to the parameters MU, SIGMA, and the mixture 
%   parameters W. MU is a D by K matrix, SIGMA is a D by D by K matrix,
%   and W is a K by 1 matrix.
%
%   If W is not provided or is equal to [], then the mixture parameters
%   are assumed to be equal in magnitude.
%
%   M
%    .Alpha
%    .Mu
%    .Sigma

% Scott Gaffney   10 June 1998
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'gmix_sim';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


% Handle preprocessing 
if (nargin < 2)
    error('error: you must supply at least 2 parameters.');
end

[D,K] = size(M.Mu);
if (~isfield(M,'Alpha') | isempty(M.Alpha))
  M.Alpha = ones(K,1)./ K;
elseif (sum(M.Alpha) ~= 1)
  error('error: the mixture weights do not add up to 1.');
end

X = zeros(N,D);
C = zeros(N,1);
Ci = zeros(K+1,1);
MLMU = zeros(D,K);
MLSIGMA = zeros(D,D,K);

% Stratify the data sets to match the given priors
C = pmfrnd(M.Alpha,N);
C = sort(C);  C=C(:);        % sort them just because
Ci = cumsum([1 hist(C,K)]);  % find class offsets

% Generate the data
for k = 1:K
  X(Ci(k):Ci(k+1)-1,:) = mvnrnd(M.Mu(:,k), M.Sigma(:,:,k), Ci(k+1)-Ci(k));
  [MLMU(:,k), MLSIGMA(:,:,k)] = normmle(X(Ci(k):Ci(k+1)-1,:));
end

M.TrueC = C;