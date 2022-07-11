function ah = SetLabelPos(ah,pos,ind,label)
%SetLabelPos  Sets the position of the axis labels appropriately
%   ah = SetLabelPos(ah,[pos],[ind],['label'])
%
%   This sets the position, in axis coordinates, of the specified
%   label.
%
% pos: offset from current value
% ind: 1 - x-pos, 2 - y-pos, 3 - z-pos
% label: 'xlabel' or 'ylabel' or 'zlabel'
%
% Defaults:
%  pos = [.25];
%  ind = 1;
%  label = 'ylabel';
%
% Example:
%  You can set the y-coordinate of the xlabel to 1 unit off set from
%  it current position as follows.
%
%  SetLabelPos(gca,-1,2,'xlabel');

% Scott Gaffney   29 January 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'SetLabelPos';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

pos = cexist('pos',.25);
ind = cexist('ind',1);
label = cexist('label','ylabel');

hnd = get(ah,label);
hpos = get(hnd,'Position');
hpos(ind) = hpos(ind) + pos;
set(hnd,'Position',hpos);
