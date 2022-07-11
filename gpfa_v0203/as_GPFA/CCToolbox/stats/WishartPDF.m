function y = WishartPDF(X,df,Sigma)
%WishartPDF  Calculate wishart density values at X.
%   WishartPDF(X,DF,Sigma) is a vector of density values calculated at X from
%   the Wishart density with parameters DF and Sigma, where DF
%   is the degress of freedom and Sigma is the scaling matrix of size
%   P-by-P. Note that DF must be greater than P-1.
%

% Scott Gaffney   14 January 2002
% DataLab @UCI
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'WishartPDF';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

[p1,p2] = size(Sigma);
if (p1~=p2)
  errorbox('Argument Error: Sigma must be square.',PROGNAME);
end
P = p1;

if (df <= P-1)
  errorbox('The degrees of freedom must be greater than P-1.',PROGNAME);
  error([PROGNAME,': The degrees of freedom must be greater than P-1.']);
end

%% Normalizing constant (its the same for the Inverse Wishart as well)
C = 2^(df*P/2)* pi^(P*(P-1)/4)* det(Sigma)^(-df/2);
for p=1:P
  C = C* gamma( (df+1-p)/2 );
end

y = (1/C)* det(X)^((df-P-1)/2)* exp(-trace(inv(Sigma)*X)/2);
