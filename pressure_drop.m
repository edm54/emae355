function [pressure_loss] = pressure_drop(m_dot, pressure1)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
    total_length = 3200;
    temp_range = [15+273 129+273];
    len = [0 3200];
    p6 = pressure1;
    L1 = 3.2 * 1000; % meters
    gravity = 9.81;
    diameter = 0.06096; % meter
    roughness = .05 ;% mm, from https://neutrium.net/fluid_flow/absolute-roughness/
    % new = .02 -.05 mm, slightly corroded = .05 - .15, moderate corrosion =
    % .15-1mm
    temp = (15 + 129)/2 + 273;
    e_d_ratio = roughness/1000 * diameter;
    area = pi * (diameter/2)^2 ;
    current_pressure = p6;
    i = 1;
    dynamic_v = refpropm('V','T',temp,'P',current_pressure/1e3, 'CO2');
    rho1 = refpropm('D','T',temp,'P',current_pressure/1e3, 'CO2');
    pressure_loss = 0;
    delta_l = 100;
    for height = 0 : delta_l : len(2)
        temp = interp1(len, temp_range, height, 'linear');
        dynamic_v = refpropm('V','T',temp,'P',current_pressure/1e3, 'CO2');
        rho(i) = refpropm('D','T',temp,'P',current_pressure/1e3, 'CO2');
        %rho(i,j) = 400
        %rho(i) = rho1;
        
        velo = m_dot/(rho(i) * area);
        reynolds = rho(i) * velo * diameter/dynamic_v;

        f1 = -1.8 * log10((6.9 / reynolds) + (e_d_ratio/ 3.7)^1.);
        f(i) = (1/f1)^2;

        head_loss(i) = (f(i) * delta_l * velo^2) / (diameter * 2  * gravity);
        pressure_drop(i) = rho(i) * gravity * head_loss(i); %pa
        pressure_loss = pressure_loss + pressure_drop(i);
        current_pressure = current_pressure + pressure_drop(i);
        i = i + 1;
    end
%      disp(current_pressure)
%      disp(pressure_loss)
end

