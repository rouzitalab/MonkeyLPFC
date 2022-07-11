function ah = SetAxesPos(ah,ind,value)
%SetAxesPos  Sets axes position values.
%   ah = SetAxesPos(ah,[ind],[value])
%
%   This function sets position(ind(i)) of the axes in ah to 
%   value(ind(i)). If value is a scalar, then 
%   it is used for each index in ind.
%
%   Defaults
%   -----------------------------
%   ind    - [1 2 3 4]
%   value  - [.25 .25 .65 .65]
%

% Scott Gaffney   29 January 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'SetAxesPos';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

ind = cexist('ind',[1 2 3 4]);
value = cexist('value',[.25 .25 .65 .65]);

value = value(:)';
if (length(value)==1)
  value = value*ones(1,length(ind));
end

for i=1:length(ah)
  pos = get(ah(i),'Position');
  pos(ind) = value;
  set(ah(i),'Position',pos);
end