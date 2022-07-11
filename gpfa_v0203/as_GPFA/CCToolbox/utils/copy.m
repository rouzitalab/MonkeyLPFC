function B = copy(A,Key)
%COPY  Copies a vector and repeats its elements according to a key or plan.
%   COPY(A,KEY) returns a vector of length sum(KEY) with the same
%   elements as A, but each element in A is repeated according to its
%   corresponding repetition number in KEY.
%
%   A can be a matrix also, in which case the dimensions of Key will
%   reference the first dimension of A.

% Scott Gaffney   14 February 2002
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'copy';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

if (prod(size(A))==length(A))
  A = A(:);
end
Key = Key(:)';
indx = cumsum([1 Key]);
B = zeros(1,sum(Key));
B(indx(1:end-1)) = 1;
trash = cumsum(B);
B = A(trash,:);
%clear trash A Key;


