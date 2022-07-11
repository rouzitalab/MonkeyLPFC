function [mu, labels] = kmeans(X, K)
%KMEANS  Perform K-means clustering on multivariate data.
%     [MU,LABELS] = KMEANS(X,K) calculates K locations for the
%     clusters that are found in the multivariate data X. The locations
%     are returned in MU, while the cluster that each data point
%     in X is assigned to is returned in LABELS.
%
%     X contains N rows of data points, each of dimension D. MU is
%     a D by K matrix, and LABELS is a N by 1 vector.

% Scott J. Gaffney	10 June 1998 
% Department of Information and Computer Science
% University of California, Irvine
%
% SJG:  26 January 1999
% Changed the handling of empty clusters to take into account
% the most likely cluster. Also changed the way that seeds (means)
% were selected for initialization.


PROGNAME = 'kmeans';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


% Preprocessing
if (K < 1)
  error('kmeans: K must be at least 1.');
end
if (isempty(X)) return, end;
empty = zeros(K,1);  % used to detect empty clusters
interval = [min(X); max(X)]';
[N,D] = size(X);
mu = zeros(D,K);
labels = zeros(N,1);

% Seed the means (or sow the fields)
pickPoints = 0;
if (pickPoints)
  distX = sum((reshape(repmat(X,1,N)',D,N*N) - repmat(X',1,N)) .^2);
  distX = reshape(distX,N,N);
  [maxX,colMax] = max(distX, [], 2);  % pick the max col for each row
  [maxX,rowMax] = max(maxX);          % pick the max row from the above cols
  mu(:,1) = X(rowMax,:)';
  mu(:,2) = X(colMax(rowMax),:)';

  selected = [rowMax, colMax(rowMax)];
  for j=3:K
    sumDist = sum(distX(selected,:));
    sumDist(selected) = 0;
    [maxX,I] = max(sumDist);
    mu(:,j) = X(I,:)';
    selected = [selected I];
  end
  clear distX;  % free this memory, we don't need it anymore.
end

pickRandom = 1;
if (pickRandom)
  perm = randperm(N);
  mu = X(perm(1:K),:)';
else     % pick equal "spaced" centers
  firstK = (interval(:,2) - interval(:,1)) / (K+1) + interval(:,1);
  for i = 1:K
    mu(:,i) = firstK * i; 
  end
end

% Make the X matrix used in computing metrics, and initialize labels.
X_NxDxK = repmat(X,[1 1 K]); 
new_labels = labels;
new_labels(1) = 1;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform K-means clustering
%
num = 0;
while (any(labels ~= new_labels))
  
  num = num + 1;
  if (rem(num,50)==0)
    fprintf('\nkmeans iteration count: %i\n',num); 
    %keyboard;
  end
  % Assign points to their nearest centers
  labels = new_labels;
  mu_NxDxK = repmat(reshape(mu,[1 D K]),[N 1 1]);
  metric_NxK = squeeze(sum((X_NxDxK - mu_NxDxK) .^2,2));
  [trash new_labels] = min(metric_NxK,[],2);  % assign new labels
  
  % Calculate the effect of the labeling on the K centers.
  for i = 1:K 
    j = find(new_labels==i);
    if (~isempty(j)),
      mu(:,i) = mean(X(j,:))';
    else
      empty(i) = 1; % We have an empty cluster (Should we decrement K?)
    end
  end
  
  % This code deals with empty clusters.
  if (any(empty))
    disp('kmeans: an empty cluster was detected.');
    
    % Find "most likely" cluster
    numlabels = hist(new_labels,1:K);
    [maxNum, maxLabel] = max(numlabels);
    indx = find(new_labels==maxLabel);
    
    % Select points from this cluster randomly
    badK = find(empty==1)';
    rperm = randperm(length(indx));
    chosen = indx(rperm(1:length(badK)));
    mu(:,badK) = X(chosen,:)';
    new_labels(chosen) = badK;
    empty(badK) = 0;
  end
end

