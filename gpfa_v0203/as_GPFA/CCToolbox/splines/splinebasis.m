function [S,S_uniq] = splinebasis(type,knots,d,X)
%SPLINEBASIS  Calculate specified spline basis functions
%
%   [S,S_uniq] = SPLINEBASIS(TYPE,D,X)
%
%     TYPE
%       'bspline'   : use B-spline basis (default)
%       'power'     : use truncated power series basis
%     D             : degree of spline (e.g., 3 for cubic)
%     X             : values at which to evalute the basis functions
%     S             : spline basis matrix
%     S_uniq        : unique rows of S (usually much smaller than S)

% Scott J Gaffney   23 January 2003
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'splinebasis';
if (~nargin)
  help(PROGNAME);
  return;
end

[uniq,trash,j] = unique(X);
switch (type)
  case 'bspline'
    S_uniq = bsplinebasis(knots,d+1,uniq);
  case 'power'
    S_uniq = powersplinebasis(knots,d,uniq);
  otherwise
    error([PROGNAME,': unrecognized basis type']);
end

S = S_uniq(j,:);

