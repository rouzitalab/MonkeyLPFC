function ah = SpreadAxes(ah,shift,offset)
%SpreadAxes  Spread the axes in a figure for better separation.
%   ah = SpreadAxes(ah,[shift],[offset])
%
%   This function spreads out (figure lengthwise) the axes that are on a single figure
%   so that there is more room between them. shift gives the amount that
%   each axes is spread out by. offset is a common length for which all
%   the axes are moved by as a group (figure lengthwise).
%
%   Defaults
%   --------------------------
%   shift    - 0.04
%   offset   - -0.05
%

% Scott Gaffney   29 January 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'SpreadAxes';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

shift = cexist('shift',0.04);
offset = cexist('offset',-0.05);

numaxes = length(ah);
for i=1:length(ah)
  AddAxesPos(ah(i),1,offset+shift*i);
end
