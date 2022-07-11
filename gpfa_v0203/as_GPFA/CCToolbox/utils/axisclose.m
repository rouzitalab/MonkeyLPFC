function ah = axisclose(ah,delta)
%AxisClose  Set axis limits to tight and then loosen just a bit.
%   ah = AxisClose(ah,[delta]), delta can also be delta = [delta(1) delta(2)]
%   where delta(1) is for the x-axis and delta(2) is for the y-axis.
%
%   This function calls AXIS AUTO and AXIS TIGHT and then adds delta to the ends
%   of each axis limit vector.
%
%   Defaults
%   ----------------------------
%   delta   - 1
%

% Scott Gaffney   29 January 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'axisclose';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

delta = cexist('delta',1);
if (prod(size(delta))==1)
  delta(2) = delta;
end

for i=1:length(ah)
  %axis(ah(i),'auto');   % make sure that AXIS TIGHT works correctly
  set(ah(i),'xlimmode','auto');
  set(ah(i),'ylimmode','auto');
  axis(ah(i),'tight');
  xlim = get(ah(i),'xlim');
  ylim = get(ah(i),'ylim');
  xlim(1) = xlim(1) - floor(delta(1)/2);
  xlim(2) = xlim(2) + floor(delta(1)/2);
  ylim(1) = ylim(1) - floor(delta(2)/2);
  ylim(2) = ylim(2) + floor(delta(2)/2);
  set(ah(i),'xlim',xlim);
  set(ah(i),'ylim',ylim);
end