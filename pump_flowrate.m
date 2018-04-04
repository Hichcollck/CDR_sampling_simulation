function [Qp] = pump_flowrate(p)

% input p in mbar

load('pumpflowrate'); % Flowrate profile from NMP 850.1.2 KNDC B datasheet
p_data = pumpflowrate(:,1);
Qp_data = pumpflowrate(:,2);

% Find the interval between the data where input z fits
j=1; % Initialization
while p>p_data(j)
    j=j+1;
end
if j==1
    Qp=0;
else
    % lineal regression
    Qp = Qp_data(j-1)-(p_data(j-1)-p)*(Qp_data(j-1)-Qp_data(j))/(p_data(j-1)-p_data(j)); % [L/min]
end

end

