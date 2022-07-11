function S = powersplinebasis(knots,d,X)
%POWERSPLINEBASIS  Calculate the truncated power series basis functions
%   S = POWERSPLINEBASIS(KNOTS,D,X)
%     KNOTS - sequence of knots
%     D     - degree of spline (e.g., 3 for cubic)
%     X     - values at which to evalute the basis functions

% Scott J Gaffney   23 January 2003
% Department of Information and Computer Science
% University of California, Irvine.

PROGNAME = 'powersplinebasis';
if (~nargin)
  help(PROGNAME);
  return;
end

% Check for knot sequence problems
if (~isempty(find(diff(knots)<0)))
  S=[];
  errorbox('Argument Error: KNOTS must be nondecreasing.',PROGNAME);
  return;
end

X = X(:);
knots = knots(:)';
m = length(knots);
n = length(X);

% generate monomial basis first
S = regmat(X,d);

% now we can add in the truncated power basis
repX = X*ones(1,m);
repKnots = ones(n,1)*knots;
pb = (repX-repKnots);
trunc = find(pb<0);
if (~isempty(trunc))
  pb(trunc) = 0;
end
S = [S pb.^d];
  
  
  
  