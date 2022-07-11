function [sorted,map] = sortcls(C,trueC)
%SORTCLS  Sort classified data without destroying classification.
%   [SORTC,MAP] = SORTCLS(C,[TrueC])

% Scott Gaffney   16 February 1999
% DataLab@UCI
% Department of Information and Computer Science
% University of California, Irvine, USA.
%
% Changes
% -------------------------------
% 5 October 2001 (Scott J Gaffney)
%   Added a warning dialog box when more than one class
%   of the true labels maps to the same class of the learned
%   labels.
%
% 12 October 2001 (Scott J Gaffney)
%   Added help/usage message on error.

PROGNAME = 'sortcls';
if (~nargin)
  help(PROGNAME);
  return;
end


if (exist('trueC','var') ~= 1)
  trueC = [];
else
  trueC = trueC(:);
end

i=1;
left = C;
map=[];
sorted = zeros(size(C));

if (isempty(trueC))
  while (~isempty(left))
    map(left(1)) = i;  % keep track of class mapping
    indx = find(C==left(1));
    sorted(indx) = i;
    left(find(left==left(1))) = [];
    i=i+1;
  end
else
  numTrueC = length(unique(trueC));
  for i=1:numTrueC
    counts(i,:) = hist(C(find(trueC==i)),1:numTrueC);
  end

  uniqC = unique(trueC);
  for i=1:numTrueC
    [trash, maxindex] = max(counts(i,:),[],2);
    map(i) = maxindex;
    sorted(find(C==maxindex) ) = i;
    counts(:,maxindex) = -1;
  end
end
  
