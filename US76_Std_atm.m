function [p,T,d] = US76_Std_atm(z)
% 1976 US Standard atmosphere (http://www.pdas.com/atmos.html)
% It contains data from 0km to 86km
load('USatm');
j=1;
while USstandardatm(j,1)<=40 % obtain the data up to 40 km
    z_US76atm(j) = USstandardatm(j,1);       % [km]
    p_US76atm(j) = USstandardatm(j,6)/100;   % Transformation from Pa to mbar
    T_US76atm(j) = USstandardatm(j,5);       % [K]
    d_US76atm(j) = USstandardatm(j,7);       % [kg/m^3]
    j=j+1;
end

% Find the interval between the data where input z fits
j=1; % Initialization
while z/1000>z_US76atm(j)
    j=j+1;
end
% lineal regression
p = p_US76atm(j-1)-(z_US76atm(j-1)-z/1000)*(p_US76atm(j-1)-p_US76atm(j))/(z_US76atm(j-1)-z_US76atm(j)); % [mbar]
T = T_US76atm(j-1)-(z_US76atm(j-1)-z/1000)*(T_US76atm(j-1)-T_US76atm(j))/(z_US76atm(j-1)-z_US76atm(j)); % [K]
d = d_US76atm(j-1)-(z_US76atm(j-1)-z/1000)*(d_US76atm(j-1)-d_US76atm(j))/(z_US76atm(j-1)-z_US76atm(j)); % [kg/m^3]

end

