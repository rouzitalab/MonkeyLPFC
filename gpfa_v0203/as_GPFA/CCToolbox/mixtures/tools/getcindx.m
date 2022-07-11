function a = getcindx(x,C)
%GETCINDX  Get the appropriate parameter index based on classification
%   GETCINDX(x,C)
%

% Scott J Gaffney   19 September 2003
% School of Computer Science
% University of California, Irvine

PROGNAME = 'getcindx';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

[n,K,D] = size(x);
for d=1:D
  indx = sub2ind([n,K,D],(1:n)',C,d(ones(n,1)));
  a(:,d) = x(indx);
end
