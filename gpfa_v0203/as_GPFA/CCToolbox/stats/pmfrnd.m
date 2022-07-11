function samp = pmfrnd(probs,N)
%PMFRND  Sample a probability mass function.
%   SAMP = PMFRND(PROBS,N) is a vector of length N containing samples
%   from the discrete probability mass function defined by
%   the vector of probabilities PROBS.

% Scott Gaffney   14 May 1998
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'pmfrnd';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

for i = 2:length(probs)
    probs(i) = probs(i-1)+probs(i);
end

samp = zeros(N,1)';
r = rand(N,1);
for i = 1:N
    for j = 1:length(probs)
	if (r(i) <= probs(j))
	    samp(i) = j;
	    break;
	end
    end
end

