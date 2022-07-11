function val = cexist(VarName,DefValue)
%CEXIST  Set default values of arguments.
%   Value = CEXIST(VarName,DefValue)

% Scott Gaffney   26 March 2002
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'cexist';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

if (evalin('caller',['exist(''',VarName,''')']) ~=1 | ...
      evalin('caller',['isempty(',VarName,')']))
  val = DefValue;
else
  val = evalin('caller',VarName);
end