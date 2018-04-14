function [z, v_z, t_cut] = SSC_algorithm(m_g,z_0,z_end)

global g

% t_cut = time elapsed since cut-off [s]
t = 1; % timestep for each iteration [s]
% v_z = vertical velocity positive up [m/s]
v_z0 = 0; % initial velocity in up direction [m/s]

C_d = 0.79; % drag coefficient for the parachute [0.75 - 0.79]
A = 80; % Surface area of parachute [m^2]

i=1;
t_cut(1)=0;
z_i=z_0;
while z_i>z_end

    [p(i), T(i), d(i)] = US76_Std_atm(z_i); % 1976 US Standard atmosphere
    v_t(i) = -((2*m_g*g)/(C_d*d(i)*A))^(1/2); % terminal velocity [m/s]

    K1 = -g*t_cut/v_t(i); % Phase coefficient
    % Phase I: canopy is considered to be uninflated (K<3)
    % Phase II: canopy is considered to be fully inflated (K>=3). The velocity
    % in this phase is close to the terminal velocity.

    if t_cut(i) < (3*v_t(i)/g) % Phase I
        v_z(i) = v_t(i)*tanh(K1);
        z_i = z_0-((v_t(i)^2)/g)*log(cosh(K1));
    else % Phase II
        K2 = v_z0/v_t(i);
        if K2<1
            v_z(i) = v_t(i)*tanh(1/(K2)-g*t/v_t(i));
            z_i = z_0 + ((v_t(i)^2)/g)*(log(cosh(1/tanh(K2)))-log(cosh(1/tanh(K2)-g*t/v_t(i))));
        elseif K2>1
            v_z(i) = v_t(i)*coth(1/coth(K2)-g*t/v_t(i));
            z_i = z_0+((v_t(i)^2)/g)*(log(sinh(1/coth(K2)))-log(sinh(1/coth(K2)-g*t/v_t(i))));
        else % K=1
            v_z(i) = v_t(i);
            z_i = z_0+v_t(i)*t;
        end
    end
    i=i+1;
    z(i)=z_i;
    t_cut(i) = t_cut(i-1)+t;
end

end

