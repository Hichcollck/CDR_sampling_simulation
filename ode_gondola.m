function [dX] = ode_gondola(t,X)
% In order to solve a second-order ODE with ode45 function
% is necessary to rewrite the equation to a system of first
% order ODEs using a change of variables

global g

% Initial conditions
Z1 = X(1);   % position (z)
Z2 = X(2);   % velocity (dz/dt)

%% Inputs
m_g = 250;      % [kg] Mass (gondola + parachute + all extra equipment)
A_g = 1.3456;   % [m^2] Bottom area (characteristic value for aerodynamic forces)
Cd_g = 1.05;    % Drag coefficient (Typical for a cube shape)
Cd_p = 1.75;    % Drag coefficient (Typical)
A_p = 80;       % [m^2] Area
C = (Cd_p*A_p+Cd_g*A_g)/m_g; % Constant

% 1976 US Standard atmosphere (http://www.pdas.com/atmos.html)
[z_USatm,p_USatm,T_USatm,d_USatm] = US76_atm(35); % obtain the data up to 35km

if Z1>1000 % "Stop" the simulation down to 1km of altitude
    d = density(Z1,z_USatm,d_USatm);
else
    d=0;
end

%% ODE
% The velocity of the satellite is given by the derivate of Z, the
% following equation is integrated from the second derivate of R
dZ2 = C*(Z2^2)*d-g;
% Then, the position vector is found by integrating the velocity obtained
dZ1 = Z2;

dX=[dZ1;dZ2];
