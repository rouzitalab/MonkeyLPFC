function ah = AddAxesPos(ah,ind,value)
%AddAxesPos  Add values to current axes positions.
%   ah = AddAxesPos(ah,ind,value)
%
%   This function sets position(ind(i)) of the axes in ah to 
%   value(ind(i)) + position(ind(i)). If value is a scalar, then 
%   it is used for each index in ind.
%

% Scott Gaffney   29 January 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'AddAxesPos';
if (nargin ~= 3)
  try; help(PROGNAME); catch; end
  return;
end

if (length(value)==1)
  value = value*ones(1,length(ind));
end

for i=1:length(ah)
  pos = get(ah(i),'Position');
  pos(ind) = pos(ind) + value;
  set(ah(i),'Position',pos);
end