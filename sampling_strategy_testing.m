%%-------------------------------------------------------------------------
%           Alternative Aircore Atmospheric Trace Gas Sampling 
%              REXUS/BEXUS 26/27 Programme (SNSB, DLR, ESA)                         
%%-------------------------------------------------------------------------
%     Air Sampling Simulator for Balloon Vertical Flight (in Artic Zone)
%             "BEXUS Air Sampling Simulator (BASS - BEXUS ASS)"
%%-------------------------------------------------------------------------           
%                            TEAM TUBULAR  
%                      Author: Jordi Coll Ortega
%                   Lulea University of Technology
%%-------------------------------------------------------------------------
%                         Date: 6th September 2018
%%-------------------------------------------------------------------------

clear all 
close all
clc

global R_spc g d_sealvl

%% Constants
R_spc = 287.058;    % [J·kg^-1·K^-1] Specific ideal gas constant for dry air
g = 9.81;           % [m/s^2] gravity
d_sealvl = 1.225;   % [kg/m^3] Air density at sea level

%% Inputs
%%----------------------------
% Flight characteristics
z_ini = 500;           % [m] Initial height to start flight simulation (Esrange)
z_float = 26000;       % [m] Float altitude
t_float = 2.5*3600;    % [sec] Floatting time
z_cutoff = z_float;    % [m] Cut-off altitude
z_end = 1000;          % [m] Final height to stop flight simulation
[p_float, T_float, d_float] = US76_Std_atm(z_float); % 1976 US Standard atmosphere

% Sampling equipment
V_min = 0.4;         % [L] (180mL) Minimum volume to collect with the bags at sea lvl
bag_volume = 3;      % [L] Physical volume of the chosen bags
p_bag_max = 137.895; % [mbar] Maximum pressure stand by the bags (2psi)
flushing_time = 60;  % [sec] Time required to flush the system during flight before one sampling

% Sampling strategy
num_bags = 6;        % Number of bags used in the experiment
bag = [1, "asc", 16000, 3;
       2, "asc", 21000, 3;
       3, "desc", 22000, 3;
       4, "desc", 18000, 3;
       5, "desc", 15000, 3;
       6, "desc", 12000, 3]; 
       % Columns: bag number, Phase (asc/desc), Opening altitude [m], Sample volume [L]

% Gondola
m_g = 250;           % [kg] Mass in DESCENT (gondola + parachute + all extra equipment)
A_g = 1.3456;        % [m^2] Bottom area (characteristic value for aerodynamic forces)
Cd_g = 1.05;         % Drag coefficient (Typical for a cube shape)
% Parachute
Cd_p = 1.75;         % Drag coefficient (Typical for parachutes)
A_p = 80;            % [m^2] Area
%%----------------------------

C = (Cd_p*A_p+Cd_g*A_g)/m_g;     % Dynamic const1ant
m_air = (V_min/1000)*d_sealvl;   % [kg] Minimum amount of air to be collected by the bags

%% Flight Profile 
% Initialization of dynamic vectors
z(1) = z_ini;
z_i = z(1);     % auxiliary copy of z variable
dt(1)=0;
time(1)=0;
v_th(1)=5;      % Initial vertical velocity
a_th(1)=0;      % Initial hypothetical acceleration
 
i=1;            % auxiliary value to create the vectors for each variable according to the altitude.
% Ascent phase
while z_i<z_float
    [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
    if i>1
    z(i)=z_i;
    v_th(i)=5; % [m/s] aproximately constant
    a_th(i)=0;
    dt(i)=(z(i)-z(i-1))/v_th(i);
    time(i)=time(i-1)+dt(i);
    end
    i=i+1;
    z_i=z(i-1)+2;    % (Meshing) Analysis each 2 meters.
end
% Float phase (assumed constant altitude)
z(i)=z_i;
time(i)=time(i-1)+t_float;
v_th(i)=0; 
a_th(i)=-10*g;
i=i+1;
% Descent phase
while z_i>z_end
    [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
    z(i)=z_i;
    [v_th(i),a_th(i),time(i),dt(i)]=free_fall_dynamics(i, z, v_th, a_th, time, d(i), C);
    i=i+1;
    z_i=z(i-1)-2;    % (Meshing) Analysis each 2 meters.
end


%% Bags sampling

% p_total(i)=p(i)+(0.5*d(i)*(v_th(i))^2)/100; % [mbar] Dynamic pressure
% [Qp(i)]=pump_flowrate(p_total(i)); % [L/min]
% V_sample(1)=0;
% V_total(1)=0;
% m_sample(1)=0;
% m_total(1)=0;
% total_bag_V(1)=0;
%  
 
for j=1:num_bags
i=1;
mass_sample=0;      % mass counter
volume_sample=0;    % volume counter
dp_sample=0;        % pressure difference counter

%%------------------- Ascent phase algorithm ------------------------------
if bag(j,2)=="asc"
z_i=z(1);
    % search opening altitude
    while bag(j,3)>z_i
       [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
       V_air(i) = m_air/d(i)*1000; % [L]
       bag_z(j,i)=NaN;
       bag_t(j,i)=NaN;
       flush_z(j,i)=NaN;
       flush_t(j,i)=NaN;
       bag_V(j,i)=0; % vacuum in bags assumed before sampling
       min_sample_V(i)= V_min;
       p_ext=p(i);
       p_int=p(i)*V_air(i)/bag_V(j,i);
       bag_dp(j,i)= p_int-p_ext;
       max_bag_p(i)= p_bag_max; 
       i=i+1;
       z_i=z(i);
    end
    % rewrite flush variables to introduce the flushing time gap before the sample starts
    [flush_z(j,i), flush_t(j,i)] = flushing(i, z, time, flush_z, flush_t, flushing_time);
    
    % It is considered that the bags will be filled when it has been
    % sampled the minimum mass required to be able to perform the analysis
    % or when the 80% of the volume of the bag is sampled (vendor recommendation).
    % while (bag not full yet) = sampling
    while (mass_sample<m_air)&&(volume_sample<(0.8*bag_volume)
       [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
       V_air(i) = m_air/d(i)*1000; % [L]
       Qp(i)=2; % [L/min] constant value 
%        p_total(i)=p(i)+(0.5*d(i)*(v_th(i))^2)/100; % [mbar] Dynamic pressure
%        [Qp(i)]=pump_flowrate(p_total(i)); % [L/min] flowrate profile depending on pump efficiency  
       bag_z(j,i)=z(i);
       bag_t(j,i)=time(i);
       bag_V(j,i)=Qp(i)*dt(i)/60; % [L];
%       total_bag_V(i)= total_bag_V(i-1)+bag_V(j,i);
       min_sample_V(i)= V_min;
       bag_p(j,i)= ;
       max_bag_p(i)= p_bag_max;
    end
    % from closing bag to end of flight    
    while time_i<time(lenght(time))
       [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
       V_air(i) = m_air/d(i)*1000; % [L]
       bag_z(j,i)=NaN;
       bag_t(j,i)=NaN;
       bag_V(j,i)= ;
       min_sample_V(i)= V_min;
       p_ext=p(i);
       p_int=p(i)*V_air(i)/bag_V(j,i);
       bag_dp(j,i)= p_int-p_ext;
       max_bag_p(i)= p_bag_max; 
       i=i+1;
       time_i=time(i);
    end    
    
%%----------------- Descent phase algorithm -------------------------------
elseif bag(j,2)=="desc"
z_i=z(1);
    % from start flight to end floating phase
    while 
        
    end
    % search of opening altitude
    while
        
    end
    % from closing bag to end flight
    while
        
    end

end

%% Total Volume Calculation
for i=1:length(time)
    total_bag_V(i) = bag_V(1,i)+bag_V(2,i)+bag_V(3,i)+bag_V(4,i)+bag_V(5,i)+bag_V(6,i);
end

%% Output plotting

figure('Name','altitude vs. time')
hold on
plot(time/3600, z/1000,'HandleVisibility','off')
scatter(bag_t(1,:)/3600, bag_z(1,:)/1000,'r','filled') % change color bags
scatter(bag_t(2,:)/3600, bag_z(2,:)/1000,'r','filled')
scatter(bag_t(3,:)/3600, bag_z(3,:)/1000,'r','filled')
scatter(bag_t(4,:)/3600, bag_z(4,:)/1000,'r','filled')
scatter(bag_t(5,:)/3600, bag_z(5,:)/1000,'r','filled')
scatter(bag_t(6,:)/3600, bag_z(6,:)/1000,'r','filled')
scatter(flush_t/3600, flush_z/1000,'y','filled')
xlabel('Time [hours]')
ylabel('Altitude [km]')
legend('Bag 1', 'Bag 2', 'Bag 3', 'Bag 4', 'Bag 5', 'Bag 6', 'Flushing')
hold off

% figure('Name','velocity vs altitude')
% plot(z/1000, v_th)
% xlabel('Altitude [km]')
% ylabel('Velocity [m/s]')
 
% figure('Name','acceleration vs altitude')
% plot(z/1000, a_th)
% xlabel('Altitude [km]')
% ylabel('Acceleration [m/s^2]')

% figure('Name','Pump performance')
% plot(Qp, z/1000)
% xlabel('Flowrate capacity [L/min]')
% ylabel('Altitude [km]')

figure('Name','Bags volume during flight')
hold on
plot(time/3600, bag_V(1,:))
plot(time/3600, bag_V(2,:))
plot(time/3600, bag_V(3,:))
plot(time/3600, bag_V(4,:))
plot(time/3600, bag_V(5,:))
plot(time/3600, bag_V(6,:))
plot(time/3600, min_sample_V)
plot(time/3600, total_bag_V)
xlabel('Time [hours]')
ylabel('Volume [L]')
legend('Bag 1', 'Bag 2', 'Bag 3', 'Bag 4', 'Bag 5', 'Bag 6', 'Total bag volume', 'Minimum sample volume')
hold off

figure('Name','Bags pressure difference during flight')
hold on
plot(time/3600, bag_p(1,:))
plot(time/3600, bag_p(2,:))
plot(time/3600, bag_p(3,:))
plot(time/3600, bag_p(4,:))
plot(time/3600, bag_p(5,:))
plot(time/3600, bag_p(6,:))
plot(time/3600, max_bag_p)
xlabel('Time [hours]')
ylabel('Pressure difference [mbar]')
legend('Bag 1', 'Bag 2', 'Bag 3', 'Bag 4', 'Bag 5', 'Bag 6', 'Bag pressure limit')
hold off

