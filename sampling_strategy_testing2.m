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
%                         Date: 10th September 2018
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
% Flight characteristics
z_ini = 500;           % [m] Initial height to start flight simulation (Esrange)
z_float = 26000;       % [m] Float altitude
t_float = 2.5*3600;    % [sec] Floatting time
z_cutoff = z_float;    % [m] Cut-off altitude
z_end = 1000;          % [m] Final height to stop flight simulation
meshing_step = 2;      % [m] (Meshing) Analysis each 2 meters.
[p_float, T_float, d_float] = US76_Std_atm(z_float); % 1976 US Standard atmosphere

% Sampling equipment
V_min = 0.4;         % [L] (180mL) Minimum volume to collect with the bags at sea lvl
bag_volume = 3;      % [L] Physical volume of the chosen bags
p_bag_max = 137.895; % [mbar] Maximum pressure stand by the bags (2psi)
flushing_time = 60;  % [sec] Time required to flush the system during flight before one sampling

% Sampling strategy
num_bags = 6;        % Number of bags used in the experiment
bag_phase = ["asc", "asc", "desc", "desc", "desc", "desc"];
bag_altitude = [18000, 22000, 24000, 18000, 14000, 10000];

% Gondola
m_g = 250;           % [kg] Mass in DESCENT (gondola + parachute + all extra equipment)
A_g = 1.3456;        % [m^2] Bottom area (characteristic value for aerodynamic forces)
Cd_g = 1.05;         % Drag coefficient (Typical for a cube shape)
% Parachute
Cd_p = 1.75;         % Drag coefficient (Typical for parachutes)
A_p = 80;            % [m^2] Area

%%-------------------------------------------------------------------------
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
%----- Ascent phase -----
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
    z_i=z(i-1)+meshing_step;
end
%----- Float phase -----
%(assumed constant altitude)
z(i)=z_i;
time(i)=time(i-1)+t_float;
v_th(i)=0; 
a_th(i)=-10*g;
i=i+1;
%----- Descent phase -----
while z_i>=z_end
    [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
    z(i)=z_i;
    [v_th(i),a_th(i),time(i),dt(i)]=free_fall_dynamics(i, z, v_th, a_th, time, d(i), C);
    i=i+1;
    z_i=z(i-1)-meshing_step;    
end

%% Bags sampling

% complete flush vectors
for i=1:length(time) 
    flush_z(i)=NaN;
    flush_t(i)=NaN;
end
 
for j=1:num_bags
i=1;
mass_sample=0;      % mass counter
volume_sample=0;    % volume counter
% dp_sample=0;        % pressure difference counter

%------------------- Ascent phase algorithm ------------------------------
if bag_phase(j)=="asc"
z_i=z(i);
    % search opening altitude
    while bag_altitude(j)>z_i
       [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
%        V_air(i) = m_air/d(i)*1000; % [L]
       bag_z(j,i)=NaN;
       bag_t(j,i)=NaN;
       
       bag_V(j,i)=0; % Assumption: 0 mbar (vacuum) in bags before sampling
       min_sample_V(i)= V_min;
       p_ext=p(i);
       p_int=0;
       bag_dp(j,i)= p_int-p_ext;
       max_bag_p(i)= p_bag_max; 
       
       i=i+1;
       z_i=z(i);
    end
    % rewrite flush variables to introduce the flushing time gap before the sample starts
    [flush_z, flush_t] = flushing(i, z, time, flush_z, flush_t, flushing_time);
    
    % It is considered that the bags will be filled when it has been
    % sampled the minimum mass required to be able to perform the analysis
    % or when the 80% of the volume of the bag is sampled (vendor recommendation).
    % while (bag not full yet) = sampling
    while (mass_sample<m_air)&&(volume_sample<(0.8*bag_volume))
       [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
       V_air(i) = m_air/d(i)*1000; % [L]
       Qp(i)=2; % [L/min] constant value 
%        p_total(i)=p(i)+(0.5*d(i)*(v_th(i))^2)/100; % [mbar] Dynamic pressure
%        [Qp(i)]=pump_flowrate(p_total(i)); % [L/min] flowrate profile depending on pump efficiency  
       bag_z(j,i)=z(i);
       bag_t(j,i)=time(i);
       bag_V_dt(j,i)=Qp(i)*dt(i)/60; % [L];
       % V2 = (p1/p2)*V1 + dV (altitude volume correction)
       volume_sample=(p(i-1)/p(i))*volume_sample+bag_V_dt(j,i); % [L]
       bag_V(j,i)=volume_sample;
       sample_m_dt(j,i)=bag_V_dt(j,i)/1000*d(i); % [kg]
       mass_sample=mass_sample+sample_m_dt(j,i);
       p_ext=p(i);
       p_int= p(i)*V_air(i)/bag_volume;
%        p_int=(mass_sample/(bag_volume/1000))*R_spc*T(i);
       bag_dp(j,i)= p_int-p_ext;
       min_sample_V(i)= V_min;
       max_bag_p(i)= p_bag_max;
        
       i=i+1;
    end
    time_i=time(i);
    % from closing bag to end of flight    
    while time_i<=time(length(time))
       [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
       V_air(i) = m_air/d(i)*1000; % [L]
       bag_z(j,i)=NaN;
       bag_t(j,i)=NaN;
       % V2 = (p1/p2)*V1 (altitude volume correction)
       bag_V(j,i)=(p(i-1)/p(i))*bag_V(j,i-1);
       p_ext=p(i);
       p_int= p(i)*V_air(i)/bag_volume;
%        p_int=(mass_sample/(bag_volume/1000))*R_spc*T(i);
       bag_dp(j,i)= p_int-p_ext;
       min_sample_V(i)= V_min;
       max_bag_p(i)= p_bag_max; 
       
       if i==length(time)
            time_i=time(i)+1;
       else
            i=i+1;
            time_i=time(i);
       end
    end    
% %----------------- Descent phase algorithm -------------------------------
elseif bag_phase(j)=="desc"
    % start again from the beginning in case no bags during ascent phase
    z_i=z(i); 
    % from start flight to floating phase
    while z_i<z_float
       [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
%        V_air(i) = m_air/d(i)*1000; % [L]
       bag_z(j,i)=NaN;
       bag_t(j,i)=NaN;
       bag_V(j,i)=0; % Assumption: 0 mbar (vacuum) in bags before sampling
       min_sample_V(i)= V_min;
       p_ext=p(i);
       p_int=0;
       bag_dp(j,i)= p_int-p_ext;
       max_bag_p(i)= p_bag_max; 
       
       i=i+1;
       z_i=z(i);
    end
    % search of opening altitude
    while z_i>bag_altitude(j)
       [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
%        V_air(i) = m_air/d(i)*1000; % [L]
       bag_z(j,i)=NaN;
       bag_t(j,i)=NaN;
       bag_V(j,i)=0; % Assumption: 0 mbar (vacuum) in bags before sampling
       min_sample_V(i)= V_min;
       p_ext=p(i);
       p_int=0;
       bag_dp(j,i)= p_int-p_ext;
       max_bag_p(i)= p_bag_max; 
       
       i=i+1;
       z_i=z(i);
    end
    % rewrite flush variables to introduce the flushing time gap before the sample starts
    [flush_z, flush_t] = flushing(i, z, time, flush_z, flush_t, flushing_time);
    
    % while (bag not full yet) = sampling
    while (mass_sample<m_air)&&(volume_sample<(0.8*bag_volume))
       [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
       V_air(i) = m_air/d(i)*1000; % [L]
       Qp(i)=2; % [L/min] constant value 
%        p_total(i)=p(i)+(0.5*d(i)*(v_th(i))^2)/100; % [mbar] Dynamic pressure
%        [Qp(i)]=pump_flowrate(p_total(i)); % [L/min] flowrate profile depending on pump efficiency  
       bag_z(j,i)=z(i);
       bag_t(j,i)=time(i);
       bag_V_dt(j,i)=Qp(i)*dt(i)/60; % [L];
       % V2 = (p1/p2)*V1 + dV (altitude volume correction)
       volume_sample=(p(i-1)/p(i))*volume_sample+bag_V_dt(j,i);
       bag_V(j,i)=volume_sample;
       sample_m_dt(j,i)=bag_V_dt(j,i)/1000*d(i);
       mass_sample=mass_sample+sample_m_dt(j,i);
       p_ext=p(i);
       p_int= p(i)*V_air(i)/bag_volume;
%        p_int=(mass_sample/(bag_volume/1000))*R_spc*T(i);
       bag_dp(j,i)= p_int-p_ext;
       min_sample_V(i)= V_min;
       max_bag_p(i)= p_bag_max;
        
       i=i+1;
    end
    time_i=time(i);
    % from closing bag to end flight
    while time_i<=time(length(time))
       [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
       V_air(i) = m_air/d(i)*1000; % [L]
       bag_z(j,i)=NaN;
       bag_t(j,i)=NaN;
       % V2 = (p1/p2)*V1 (altitude volume correction)
       bag_V(j,i)=(p(i-1)/p(i))*bag_V(j,i-1);
       p_ext=p(i);
       p_int= p(i)*V_air(i)/bag_volume;
%        p_int=(mass_sample/(bag_volume/1000))*R_spc*T(i);
       bag_dp(j,i)= p_int-p_ext;
       min_sample_V(i)= V_min;
       max_bag_p(i)= p_bag_max; 
       
       if i==length(time)
            time_i=time(i)+1;
       else
            i=i+1;
            time_i=time(i);
       end
    end
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
scatter(bag_t(1,:)/3600, bag_z(1,:)/1000,'filled') % change color bags
scatter(bag_t(2,:)/3600, bag_z(2,:)/1000,'filled')
scatter(bag_t(3,:)/3600, bag_z(3,:)/1000,'filled')
scatter(bag_t(4,:)/3600, bag_z(4,:)/1000,'filled')
scatter(bag_t(5,:)/3600, bag_z(5,:)/1000,'filled')
scatter(bag_t(6,:)/3600, bag_z(6,:)/1000,'filled')
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
plot(time/3600, bag_dp(1,:))
plot(time/3600, bag_dp(2,:))
plot(time/3600, bag_dp(3,:))
plot(time/3600, bag_dp(4,:))
plot(time/3600, bag_dp(5,:))
plot(time/3600, bag_dp(6,:))
plot(time/3600, max_bag_p)
xlabel('Time [hours]')
ylabel('Pressure difference [mbar]')
legend('Bag 1', 'Bag 2', 'Bag 3', 'Bag 4', 'Bag 5', 'Bag 6', 'Bag pressure limit')
hold off

