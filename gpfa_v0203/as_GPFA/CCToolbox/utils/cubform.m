function x = cubform(coef)
%CUBFORM  Find roots of cubic equation using the Cubic Formula
%   X = CUBFORM(coef) returns all three roots of B(4)^3 + ... + B(1). A real
%   root is always returned in X(1).

% Scott Gaffney   10 April 2003
% Department of Information and Computer Science
% University of California, Irvine

PROGNAME = 'cubform';
if (~nargin)
  try; help(PROGNAME); catch; end
  return;
end


x = zeros(3,1);
a=coef(1); b=coef(2); c=coef(3); d=coef(4);

% if delta2 <= 0, y has no extrema
delta2 = (b^2-3*a*c)/(9*a^2);

% if e2 < h2, y has 3 real solns
% if e2 > h2, y has 1 real soln
% if e2 == h2, y has 1 real repeated soln and 1 other real soln
h2 = 4*a^2*delta2^3;
e = (d + 2*b^3/(27*a^2) - b*c/(3*a));
e2 = e^2;
r = -b/(3*a);

% check for e2 == h2  (repeated soln)
if (abs(e2-h2) < eps)
  delta = sqrt(delta2);
  if (e<0)
    delta = -delta;
  end
  x(1) = delta + r;
  x(2) = x(1);
  x(3) = -2*delta + r;
  return;
end

% Solution to depressed cubic
t3 = (e + sqrt(e2 - h2))/(2*a);
s3 = t3 -e/a;  %  s3 = (-e + sqrt(e2 - h2))/(2*a)

% Get 'real' cube roots
if (t3<0)
  t = -(-t3)^(1/3);
else
  t = (t3)^(1/3);
end
if (s3<0)
  s = -(-s3)^(1/3);
else
  s = (s3)^(1/3);
end

% complex multipliers used to get complex/real solns
p = (-1 + sqrt(-3))/2;
q = conj(p);

% calculate the solns
x(1) =   s -   t + r;
x(2) = p*s - q*t + r;
x(3) = q*s - p*t + r;

% if all solns are real, then strip insignificant 0*i
if (e2 < h2)
  x = real(x);
end