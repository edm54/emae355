function [current_press] = pressure_drop_down_one(m_dot, pressure1, delta_l, temp)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    gravity = 9.81;
    diameter = 0.06096; % meter
    roughness = .05 ;% mm, from https://neutrium.net/fluid_flow/absolute-roughness/
    % new = .02 -.05 mm, slightly corroded = .05 - .15, moderate corrosion =
    % .15-1mm
    e_d_ratio = roughness/1000 * diameter;
    area = pi * (diameter/2)^2 ;
    current_pressure = pressure1;
    cp2 = pressure1;
    i = 2;
    R = 188.92; % Ideal gas Constant
    rho1 = refpropm('D','T',temp,'P',current_pressure/1e3, 'CO2');    
        
    dynamic_v1 = refpropm('V','T',temp,'P',current_pressure(i-1)/1e3, 'CO2');
    rho1(i) = refpropm('D','T',temp,'P',current_pressure(i-1)/1e3, 'CO2');

    velo1 = m_dot/(rho1(i) * area);
    reynolds(i) = rho1(i) * velo1 * diameter/dynamic_v1;

    if reynolds(i) >= 4000
        f1 = -1.8 * log10((6.9 / reynolds(i)) + (e_d_ratio/ 3.7)^1.); % turb only
        f(i) = (1/f1)^2;
    else
        f(i) = 64/reynolds(i)
    end

    head_loss(i) = (f(i) * delta_l * velo1^2) / (diameter * 2  * gravity);
    pressure_drop(i) = rho1(i) * gravity * head_loss(i); %pa
    pressure_loss = pressure_drop(i);
    current_pressure(i) = current_pressure(i-1) - pressure_drop(i);
    rho2 = refpropm('D','T',temp,'P',current_pressure(i)/1e3, 'CO2');
    velo2 = m_dot/(rho2 * area);

    velo_diff = (velo2^2 - velo1^2)/2;
    % Gravity
    %velo_diff = 0; % If we do small enough where rho does not change signifcantly(?)
    % Negative delta l since 2 is below 1
    p2 = current_pressure(i) * exp((-1/(R * temp)) * (velo_diff +  gravity * -1 * delta_l));
    gravity_gain = p2 - current_pressure(i) ;
    cp2(i) = rho1(i) * ((velo1^2 - velo2^2) * 1/2 + gravity * delta_l) - pressure_drop(i); 
    current_pressure(i) = p2;
    current_press = current_pressure(i)
end
