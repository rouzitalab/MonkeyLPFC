function trajs = fvec2cell(y)
%FVec2CELL  Convert curves in feature vector format to Cell format
%   TRAJS = FVec2CELL(Y) converts the n feature vectors in the n-by-m-by-D 
%   matrix y into a cell array of trajectories.
%

% Scott J Gaffney   06 June 2002
% Department of Information and Computer Science
% University of California, Irvine.
%


PROGNAME = 'FVec2cell';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

[n,m,D] = size(y);
for i=1:n
  trajs{i,1} = permute(y(i,:,:),[2 3 1]);
end
