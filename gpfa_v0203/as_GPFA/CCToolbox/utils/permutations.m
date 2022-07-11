function perms = permutations(x)
%permutations  Generate all permutations recursively
%   permutations(x) generates the matrix of permutations of
%   the vector x. The permutations are stored along the rows.

% Scott Gaffney   1 December 2003
% Department of Information and Computer Science
% University of California, Irvine


PROGNAME = 'permutations';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

len = length(x);
perms = zeros(prod(1:len),len);  % see !! below (we don't have to pre-allocate)
if (len==1)
  perms = x;
elseif (len==2)
  perms = [x(1) x(2); x(2) x(1)];
elseif (len==3)
  perms = [x(1) x(2) x(3);
           x(1) x(3) x(2);
           x(2) x(1) x(3);
           x(2) x(3) x(1);
           x(3) x(1) x(2);
           x(3) x(2) x(1)];
else
  n = prod(1:len-1);
  for i=1:len
    tmp = x;
    tmp(i)=[];
    p = permutations(tmp);
    %perms = [perms; ones(n,1)*x(i) p];  % !! we can build up the memory
    perms((i-1)*n+1:i*n,1) = x(i);  % store the held-out element
    perms((i-1)*n+1:i*n,2:end) = p; % store the trailing permutations
  end
end
