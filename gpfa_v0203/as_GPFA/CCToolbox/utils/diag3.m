function y = diag3(x)
%DIAG3  Get the diagonal elements for a 3-D matrix
%   DiagMatrix = DIAG3(X) returns the diagonal elements for
%   each 2-D matrix X(:,:,i) in the columns of DiagMatrix for each i .

% Scott Gaffney   04 April 2002
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'diag3';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

[r,c,n] = size(x);
if (r~=c)
  error([PROGNAME,': X must be a square matrix.']);
end
if (r==1)      % handle special vector case
  y = x(:)';
  return;
end

diagi = (1 + (r*(0:r-1)) + (0:r-1))' * ones(1,n);
diagi = diagi + ones(r,1) * (r^2*(0:n-1));
y = x(diagi);