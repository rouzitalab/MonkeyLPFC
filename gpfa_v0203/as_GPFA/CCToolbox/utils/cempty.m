function val = cempty(CurValue,DefValue)
%CEMPTY  If argument empty, return []; else return default value.
%   Value = CEMPTY(CurValue,DefValue)

% Scott Gaffney   05 April 2002
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'cempty';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

if (isempty(CurValue))
  val = [];
else
  val = DefValue;
end