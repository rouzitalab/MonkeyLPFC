function knots = selectknots(x,basis,degree,m)
%SELECTKNOTS  Select knots for spline basis.
%   S = SELECTKNOTS(x,basis,degree,multiplier)
%
%   x               : min and max limits or sample of possible x-values
%   basis
%     'bspline'     : select knots for a B-spline basis
%     'power'       : selects knots for a truncated power series spline basis
%   degree          : degree of fit between knots (e.g., 3 for cubic splines)
%   multiplier      : number multiplied by length of x-interval giving number
%                     of internal knots to place along the x-axis
%

% Scott J Gaffney   23 January 2003
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'selectknots';
if (~nargin)
  help(PROGNAME);
  return;
end


uniq = unique(x);
knots = linspace(uniq(1),uniq(end),ceil((uniq(end)-uniq(1))*m));

switch (basis)
  case 'bspline'
    knots = addendpts(knots,degree+1);
    
  case 'power'
    knots(1)  = (knots(1)+knots(2))/2;
    knots(end) = (knots(end-1)+knots(end))/2;
end
  
  
  