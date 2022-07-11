function [knots,num] = addendpts(knots,K)
%AddEndPts Add necessary endpoint knots based on order of spline
%
%   [knts,num] = AddEndPts(Knots,Order) returns the knot sequence Knots with 
%   added endpoint knots satisfying the necessary end conditions (as per
%   the specified order). The number of knots added (or subtracted) to the
%   front of Knots is returned as num.

% Scott J Gaffney   29 January 2002
% Department of Information and Computer Science
% University of California, Irvine.

if (prod(size(K))~=1),  error('The spline order must be a scalar'); end

kdelta = diff(knots);
if (~isempty(find(kdelta<0)))
  knots = sort(knots); 
  kdelta = diff(knots); 
end

changepoints = find(kdelta>0);
if (isempty(changepoints))
   error('The specified knot sequence must have more than one point')
end
 
num = K - changepoints(1);
if (length(changepoints)==1)
  interior = [];
else
  interior = [changepoints(2):changepoints(end)];
end
k1 = knots(1);
kn = knots(end);
knots = [k1(ones(1,K)) knots(interior) kn(ones(1,K))];
