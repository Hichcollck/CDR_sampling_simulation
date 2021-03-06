%%-------------------------------------------------------------------------
%           Alternative Aircore Atmospheric Trace Gas Sampling 
%              REXUS/BEXUS 26/27 Programme (SNSB, DLR, ESA)                         
%%-------------------------------------------------------------------------
%     Air Sampling Simulator for Balloon Vertical Flight (in Artic Zone)
%                    "BEXUS Air Sampling Simulator"
%%-------------------------------------------------------------------------           
%                            TEAM TUBULAR  
%                      Author: Jordi Coll Ortega
%                   Lulea University of Technology
%%-------------------------------------------------------------------------
%                        Date: 12th March 2018
%%-------------------------------------------------------------------------

clear all 
close all
clc

global R_spc R_ideal g d_sealvl

%% Constants
R_spc = 287.058;    % [J�kg^-1�K^-1] Specific ideal gas constant for dry air
R_ideal = 8314;    % Ideal gas constant
g = 9.81;           % [m/s^2] gravity
d_sealvl = 1.225;   % [kg/m^3] Air density at sea level

%% Inputs
%%----------------------------
% Flight characteristics
phase = 2;           % ascent = 1; descent = 2; 
    % If ascent phase
z_float = 26000;     % [m] Float altitude
z_ini = 18000;       % [m] Initial height to start sampling
    % If descent phase
z_cutoff = 26000;    % [m] Cut-off altitude
z_end = 12000;       % [m] Final height where the sampling stops
% Sampling equipment
V_min = 0.4;         % [L] (180mL) Minimum volume to collect with the bags at sea lvl
bag_volume = 3;      % [L] Physical volume of the chosen bags
p_bag_max = 137.895; % [mbar] Maximum pressure stand by the bags (2psi)
flushing_time = 60;  % [sec] Time required to flush the system between the sampling of two bags
% Sampling strategy
num_bags = 6;          % Number of bags used in the experiment
% z_bags = [, , , , , ]; % [m] Opening height of each bag
% Gondola
m_g = 250;           % [kg] Mass in DESCENT (gondola + parachute + all extra equipment)
A_g = 1.3456;        % [m^2] Bottom area (characteristic value for aerodynamic forces)
Cd_g = 1.05;         % Drag coefficient (Typical for a cube shape)
% Parachute
Cd_p = 1.75;         % Drag coefficient (Typical for parachutes)
A_p = 80;            % [m^2] Area
%%----------------------------

C = (Cd_p*A_p+Cd_g*A_g)/m_g;     % Dynamic constant
m_air = (V_min/1000)*d_sealvl;   % [kg] Minimum amount of air to be collected by the bags

M_air = R_ideal/R_spc ; % air molecular mass
n = m_air/M_air; % number of molecules

%% Processing
%%------------------- Ascent phase algorithm ------------------------------
if phase==1
z(1)=z_ini;       % Height vector definition
z_i=z_ini;        % auxiliary copy of z variable
i=1;              % auxiliary value to create the vectors for each variable according to the altitude.

[p_float, T_float, d_float] = US76_Std_atm(z_float); % 1976 US Standard atmosphere

% Initialize counters
m_counter=0;
bag_counter=0;
V_counter=0;

while z_i<=z_float  % Analysis only up to 1000 meters above sea level.
    [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
    V_bag(i) = m_air/d(i)*1000; % [L]
    
    Qp(i)=1; % [L/min] constant value
    if i>1
        z(i)=z_i;
        v_th(i)=5; % [m/s] aproximately constant
        a_th(i)=0;
        dt(i)=(z(i)-z(i-1))/v_th(i);
        time(i)=time(i-1)+dt(i);
        
%         p_total(i)=p(i)+(0.5*d(i)*(v_th(i))^2)/100; % [mbar] Dynamic pressure
%         [Qp(i)]=pump_flowrate(p_total(i)); % [L/min] flowrate profile depending on pump efficiency  
        
        V_sample(i)=Qp(i)*dt(i)/60; % [L]
        V_total(i)=V_total(i-1)+V_sample(i); % [L]
        m_sample(i)=(V_sample(i)/1000)*d(i); % [kg]
        m_total(i)=m_total(i-1)+m_sample(i);
    else 
        % Initialization of dynamic vectors
        dt(1)=0;
        time(1)=0;
        v_th(1)=5;       % Initial descendent velocity
        a_th(1)=0;      % Initial hypothetical acceleration
        
%         p_total(i)=p(i)+(0.5*d(i)*(v_th(i))^2)/100; % [mbar] Dynamic pressure
%         [Qp(i)]=pump_flowrate(p_total(i)); % [L/min]
        V_sample(1)=0;
        V_total(1)=0;
        m_sample(1)=0;
        m_total(1)=0; 
    end
    
    m_counter=m_counter+m_sample(i);
    V_counter=V_counter+V_sample(i);
    
    % It is considered that the bags will be filled when it has been
    % sampled the minimum mass required to be able to perform the analysis
    % or when the 80% of the volume of the bag is sampled (vendor recommendation).
    if (m_counter>m_air)&&(V_counter>=(0.8*bag_volume))
        bag_counter=bag_counter+1;
        m_counter=m_counter-m_air;
        V_counter=0;

        bag_z(i)=z_i; % save altitude when each bag is completely filled
        bag_t(i)=time(i); % save time when each bag is completely filled
    else
        bag_z(i)=NaN;
        bag_t(i)=NaN;
    end
    
    % Pressure inside the bags
    p_ext=p(i);
    p_int= p(i)*V_bag(i)/bag_volume;
  %  p_int=(m_counter*R_spc*T(i)/(bag_volume/1000))/100; % [mbar]
    dp_bag(i)=p_int-p_ext;
    dp_max(i)=p_bag_max;
    dp_float(i)=p_int-p_float;
    
    i=i+1;
    z_i=z(i-1)+2;    % (Meshing) Analysis each 2 meters.
end

%%----------------- Descent phase algorithm -------------------------------
elseif phase==2
z(1)=z_cutoff;     % Height vector definition
z_i=z_cutoff;      % auxiliary copy of z variable
i=1;               % auxiliary value to create the vectors for each variable according to the altitude.

% Initialize counters
m_counter=0;
bag_counter=0;
V_counter=0;
last_bag=0;

while z_i>=z_end  % Analysis only up to 1000 meters above sea level.
    [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
    V_bag(i) = m_air/d(i)*1000; % [L]
    
    Qp(i)=3; % [L/min] constant value
    if i>1
        z(i)=z_i;
        [v_th(i),a_th(i),time(i),dt(i)]=free_fall_dynamics(i, z, v_th, a_th, time, d(i), C);
        
%        p_total(i)=p(i)+(0.5*d(i)*(v_th(i))^2)/100; % [mbar] Dynamic pressure
%        [Qp(i)]=pump_flowrate(p(i)); % [L/min] flowrate profile depending on pump efficiency  
        if time(i)<(last_bag+flushing_time)
            Qp(i)=0; % [L/min] constant value
        end
        V_sample(i)=Qp(i)*dt(i)/60; % [L]
        V_total(i)=V_total(i-1)+V_sample(i); % [L]
        m_sample(i)=(V_sample(i)/1000)*d(i); % [kg]
        m_total(i)=m_total(i-1)+m_sample(i);
    else 
        % Initialization of dynamic vectors
        dt(1)=0;
        time(1)=0;
        v_th(1)=0;       % Initial descendent velocity
        a_th(1)=-10*g;      % Initial hypothetical acceleration
        
%         p_total(i)=p(i)+(0.5*d(i)*(v_th(i))^2)/100; % [mbar] Dynamic pressure
%         [Qp(i)]=pump_flowrate(p_total(i)); % [L/min]
        V_sample(1)=0;
        V_total(1)=0;
        m_sample(1)=0;
        m_total(1)=0; 
    end
    
    m_counter=m_counter+m_sample(i);
    V_counter=V_counter+V_sample(i);
    
    % It is considered that the bags will be filled when it has been
    % sampled the minimum mass required to be able to perform the analysis
    % or when the 80% of the volume of the bag is sampled (vendor recommendation).
    if (m_counter>m_air)&&(V_counter>=(0.8*bag_volume))
        bag_counter=bag_counter+1;
        m_counter=m_counter-m_air;
        V_counter=0;

        bag_z(i)=z_i; % save altitude when each bag is completely filled
        bag_t(i)=time(i); % save time when each bag is completely filled
        last_bag=time(i); % auxiliary var to introduce flushing time gaps
    else
        bag_z(i)=NaN;
        bag_t(i)=NaN;
    end

    % Pressure inside the bags
    p_ext= p(i);
    p_int= p(i)*V_bag(i)/bag_volume;
  %  p_int=(m_counter*R_spc*T(i)/(bag_volume/1000))/100; % [mbar]
    dp_bag(i)=p_int-p_ext;
    dp_max(i)=p_bag_max;
    
    i=i+1;
    z_i=z(i-1)-2;    % (Meshing) Analysis each 2 meters.
end
end

%%-------------------------- SSC Algorithm --------------------------------

% [z_SSC, v_SSC, time_SSC]=SSC_algorithm(m_g, z_cutoff, z_end);

%% Post-processing (Sampling Strategy)
 
% for i=1:num_bags
%     z_bags(i)
% end

%% Output plotting

figure('Name','velocity vs altitude')
plot(z/1000, v_th)
xlabel('Altitude [km]')
ylabel('Velocity [m/s]')

figure('Name','acceleration vs altitude')
plot(z/1000, a_th)
xlabel('Altitude [km]')
ylabel('Acceleration [m/s^2]')

figure('Name','altitude vs. time')
plot(time/60, z/1000)
xlabel('Time [min]')
ylabel('Altitude [km]')

% figure('Name','SSC altitude vs. time')
% plot(time_SSC/60, z_SSC/1000)
% xlabel('Time [min]')
% ylabel('Altitude [km]')

figure('Name','Pump performance')
plot(Qp, z/1000)
xlabel('Flowrate capacity [L/min]')
ylabel('Altitude [km]')

figure('Name','Bags minimum volume over altitude')
plot(V_bag, z/1000)
xlabel('Volume [L]')
ylabel('Altitude [km]')

figure('Name','Total volume collected over altitude')
hold on
plot(V_total, z/1000,'b')
scatter(V_total, bag_z/1000,'r','filled')
xlabel('Volume [L]')
ylabel('Altitude [km]')
legend('Total Volume collected','Filled bags')
hold off

figure('Name','Total volume collected over time')
hold on
plot(time/60, V_total,'b')
scatter(bag_t/60, V_total,'r','filled')
xlabel('Time [min]')
ylabel('Volume [L]')
legend('Total Volume collected','Filled bags')
hold off

figure('Name','Pressure profile inside the bags')
hold on
if(phase==1)
    plot(dp_bag, z/1000,'b')
    plot(dp_float, z/1000,'g')
else
    plot(dp_bag, z/1000,'b')
end
plot(dp_max, z/1000,'r')
xlabel('Pressure [mbar]')
ylabel('Altitude [km]')
legend('Sampling Bag Pressure','Estimated Pressure at floating','Maximum Bag Pressure')
hold off

%% Saving the workspace for "pre_graphs.m" main inputs
% save('250kg_simulation.mat')


