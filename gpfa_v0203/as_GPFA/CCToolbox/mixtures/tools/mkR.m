function R = mkR(P,alpha,beta,B)
%MKR  Build default structure governing the estimation of R.
%   R = MKR(P,[ALPHA],[BETA],[B]) returns the default structure
%   that governs the estimation of R for use with Hierarchical
%   mixtures. 
%
%   P gives the square dimensions
%   of R, that is [P,P] = size(R). ALPHA gives the degrees of
%   freedom for the Wishart distribution and BETA gives the
%   multiplicand for the diagonal of a diagonal matrix prior
%   for the Wishart distribution.  You may specify a general
%   form for the matrix parameter of the Wishart distribution
%   in B if you do not wish to use a diagonal parameter as with BETA.
%   
%   MKR(P) creates a zero prior structure where there
%   is only one covariance matrix shared among all the
%   clusters and this matrix is diagonal.

% Scott Gaffney   1 February 2002
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'mkR';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

if (exist('alpha')~=1 | isempty(alpha))
  alpha = P+1;   % no prior at all
end
if (alpha < P-1)
  errorbox('Alpha must be greater than P-1.',PROGNAME);
  error([PROGNAME,': Alpha must be greater than P-1.']);
end
if (exist('beta')~=1 | isempty(beta))
  beta = 0;    % no prior at all
end
if (exist('B')~=1 | isempty(B))
  B = diag(ones(P,1)) * beta;  % diagonal prior if beta~=0
end


R.Share = 1;
R.Diagonal = 1;
R.NoPrior = 0;

R.Wish_a = alpha;
R.Wish_iB = B;
if (~any(any(B)))
  R.NoPrior = 1;
  R.Wish_B = B;
else
  R.Wish_B = inv(B);
end