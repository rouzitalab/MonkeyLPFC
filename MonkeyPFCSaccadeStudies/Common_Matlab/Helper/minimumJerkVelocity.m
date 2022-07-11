
function [v_out] = minimumJerkVelocity(t, duration, x_start, x_stop)
%http://shadmehrlab.org/book/minimum_jerk/minimumjerk.htm

k0 = x_stop - x_start;
k1 = k0*10*( (1/duration)^3 );
k2 = -k0*15*( (1/duration)^5 );
k3 = k0*6*( (1/duration)^5 );
x_out = x_start + k1*(t.^3) + k2*(t.^4) + k3*(t.^5);
v_out = x_out;

ymin_jerk = Min_jerk([x_start x_stop 0 0 0 0], duration, 1000);

%v_out = (x_stop - x_start)*( 30*((t/duration).^2) - 60*((t/duration).^3) + 30*((t/duration).^4));
end

function output = min_jerk(xi, xf, d, t)

    output = xi + (xf - xi) .* (10 * (t/d).^3 - 15 * (t/d).^4 + 6 * (t/d).^5);

end

function ymin_jerk = Min_jerk(S,tf,fs)
%       Mehran EMADI ANDANI
% S = [PositionStart PositionFinal Velocity0 VelocityFinal Acceleration0 AccelerationFinal]
% tf = duration of the movement (sec)
% fs = sampling frequency (Hz)
% Example: 
%   tf = 0.5;
%   fs = 1000; 
%   Time = 0:1/fs:tf;
%   ymin_jerk = Min_jerk([0 1 0 0 0 0],tf,fs);
%   plot(Time,ymin_jerk);

    T = [1;1];
    for i=1:3-1
        T = [T;tf^i;tf^i];
    end 
    ymin_jerk = (S*diag(T))*ME(3,tf,fs);
end

function [yt,y] = ME(n,tf,fs)
    syms t real

    Time = 0:1/fs:tf;

    for i=0:2*n-1
        Xt(i+1,1) = t^i;
    end
    for i=1:n-1
        eval(['d' num2str(i) 'Xt = diff(Xt,i);']);
    end
    X0 = subs(Xt,0);
    for i = 1:n-1
        eval(['d' num2str(i) 'X0 = subs(d' num2str(i) 'Xt,0);']);
    end    
    Q(:,1:2) = [X0 Xt];
    for i=1:n-1,
        Q(:,2*i+1) = eval(['d' num2str(i) 'X0;']);
        Q(:,2*i+2) = eval(['d' num2str(i) 'Xt;']);
    end

    Qinv = inv(Q);

    T = [1;1];
    for i=1:n-1
        T = [T;tf^i;tf^i];
    end 
    y  = subs(Qinv,tf)*Xt./T;
    yt = subs(y,Time); % Movement Elements

end