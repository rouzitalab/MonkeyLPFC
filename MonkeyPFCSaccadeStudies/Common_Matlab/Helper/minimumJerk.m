function [output] = minimumJerk(t, duration, xi, xf, degree)
%http://shadmehrlab.org/book/minimum_jerk/minimumjerk.htm
if nargin < 5
    degree = 0;
end

t_d = t./duration;

% Default - no derivative
% output = xi + (xf - xi) .* (10*(t_d.^3) - 15*(t_d.^4) + 6*(t_d.^5));
k0 = xi;
x1 = 10*(t_d.^3);
x2 = -15*(t_d.^4);
x3 = 6*(t_d.^5);

if degree == 1
    % Velocity
    k0 = 0;
    x1 = 30*(t_d.^2);
    x2 = -60*(t_d.^3);
    x3 = 30*(t_d.^4);
elseif degree == 2
    % Acceleration
    k0 = 0;
    x1 = 60*t_d;
    x2 = -180*(t_d.^2);
    x3 = 120*(t_d.^3);
end

output = k0 + (xf - xi) .* (x1 + x2 + x3);