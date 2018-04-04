function [v_i, a_i, time_i, dt] = free_fall_dynamics(i, z, v, a, time, d, C)

global g

%% Free fall equations
sqrt = (v(i-1))^2-4*a(i-1)*(z(i-1)-z(i));
if sqrt<0
    sqrt=0;
end

dt = (3/5)*(-v(i-1)-(sqrt)^(1/2))/a(i-1); % [seconds] Time required to go from z(i-1) to z(i)
time_i = time(i-1)+dt; % [seconds] Total time from the cut-off until z(i)

v_i = v(i-1) + a(i-1)*dt; % [m/s]

%% Force equilibrium (new acceleration for next loop pass)
% Drag is divided by "Real-Life factor" equal to 3
a_i = -g+d*v_i^2*C/(2*3); % [m/s^2]

end

