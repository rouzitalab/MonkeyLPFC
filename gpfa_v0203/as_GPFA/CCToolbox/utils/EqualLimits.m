function ah = EqualLimits(ah,delta)
%EqualLimits  Set all axes to have corresponding equivalent limit vectors.
%   ah = EqualLimits(ah,[delta])
%
%   This function determines the x and y maximum limit vectors
%   across all axes in ah and then makes these two vectors common
%   across all axes.
%
%   If delta is provided, then axisclose(ah(i),delta) is called
%   for each axis first.
%

% Scott Gaffney   29 January 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'EqualLimits';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end

delta = cexist('delta',[]);
mxlim = [inf -inf];
mylim = [inf -inf];

% Calculate the maximum limit vectors
for i=1:length(ah)
  if (~isempty(delta))
    axisclose(ah(i),delta);
  end
  xlim = get(ah(i),'xlim');
  ylim = get(ah(i),'ylim');

  % Find largest xlim
  if (xlim(1) < mxlim(1))
    mxlim(1) = xlim(1);
  end
  if (xlim(2) > mxlim(2))
    mxlim(2) = xlim(2);
  end
  
  % Find largest ylim
  if (ylim(1) < mylim(1))
    mylim(1) = ylim(1);
  end
  if (ylim(2) > mylim(2))
    mylim(2) = ylim(2);
  end
end


% Now make them common
for i=1:length(ah)
  set(ah(i),'xlim',mxlim);
  set(ah(i),'ylim',mylim);
end