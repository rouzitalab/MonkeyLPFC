function [map,Acc] = MatchClusters(C,C2,K,approx)
%MatchClusters  Solve factorial class/label mapping problem when clustering.
%   [Map,Acc] = MatchClusters(C,C2,K) provides the mapping of the class labels
%   in C that best match the class labels in C2. C and C2 are class label
%   vectors, and K gives the number of unique labels that may appear in 
%   C and C2.
%
%   find(C2==i) is the best matching cluster to find(C==Map(i)).
%   In other words, if Map = [2 3 1], then the second cluster in C best
%   matches the first cluster in C2, the third cluster in C best matches
%   the second cluster in C2, and so forth.
%
%   Acc provides the classification accuracy that results from the mapping.
%
%   [Map,Acc] = MatchClusters(C,C2,K,'greedy') performs the optimization
%   using a greedy algorithm that selects each mapping separately instead
%   of jointly. Often this is equivalent to the exact factorial
%   solution. Here is an example of when they are not equivalent.
%
%   C_true = [1 1 1 1 1 2 2 2 2 2 2 2 3 3 3];
%   C      = [3 1 1 1 2 1 1 1 1 2 2 2 2 3 3];
%   [map,Acc] = MatchClusters(C,C_true,3)
%   [gmap,gAcc] = MatchClusters(C,C_true,3,'greedy')
%
%   The correct map is [1 2 3] with accuracy 53.33%. However, the
%   greedy algorithm returns the map [2 1 3] with accuracy 46.67%.
%   The reason is because C cluster 1 matches more elements
%   of C_true cluster 2 than its own cluster (C_true cluster 1).
%   This incorrect matching throws off the first step of the
%   greedy algorithm. (The problem may be more with the clustering
%   than with the matching in this case.)

% Scott Gaffney   1 December 2003
% Department of Information and Computer Science
% University of California, Irvine


PROGNAME = 'MatchClusters';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

if (exist('approx') & strcmp(approx,'greedy'))
  approx = 1;
else
  approx = 0;
end


% Exact factorial solution
if (~approx)
  pmap = permutations(1:K);
  NumPerms = size(pmap,1);
  score = zeros(1,NumPerms);
  for i=1:NumPerms
    for j=1:K
      score(i) = score(i) + sum( (C==pmap(i,j))&(C2==j) );
    end
  end
  [maxx,where] = max(score);
  map = pmap(where,:);
  Acc = maxx/length(C);

% Approximate greedy solution  
else
  for i=1:K
    for j=1:K
      score(i,j) = sum( (C==i)&(C2==j) );
    end
  end
  
  % choose mapping in greedy manner
  Acc = 0;
  for i=1:K
    maxvals = max(score,[],2);  % find max 'single' mappings
    [maxx,f] = max(maxvals);    % find overall max 'single' mapping
    match = find(score(f,:)==maxvals(f));  % find mapping for this index
    map(match) = f;  % store the mapping
    score(f,:) = -inf;  % remove this index from consideration
    score(:,match) = -inf; % remove this mapping from consideration
    Acc = Acc + maxx;
  end
  Acc = Acc/length(C);
end
