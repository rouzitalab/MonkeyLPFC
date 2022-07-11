function ah = axisequal(ah)
%AxisEqual  Set DataAspectRatio so that x,y and z axes are on the same scale.
%   ah = AxisEqual(ah)
%
%   This function sets the DataAspectRatio property to [1 1 1] so that
%   the units on all of the axes are on the same scale. If ah is not 
%   provided, then the current axes are used.
%

% Scott Gaffney   29 January 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'AxisEqual';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


set(ah,'DataAspectRatio',[1 1 1]);