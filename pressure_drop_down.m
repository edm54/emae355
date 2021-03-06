function [current_press, pressure_loss, gravity_gain, temp] = pressure_drop_down(m_dot, pressure1, t_init)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
   
    len = [0 3200];
    L1 = 3.2 * 1000; % meters
    gravity = 9.81;
    diameter = 0.06096; % meter
    roughness = .05 ;% mm, from https://neutrium.net/fluid_flow/absolute-roughness/
    % new = .02 -.05 mm, slightly corroded = .05 - .15, moderate corrosion =
    % .15-1mm
    e_d_ratio = roughness/1000 * diameter;
    area = pi * (diameter/2)^2 ;
    current_pressure = pressure1;
    cp2 = pressure1;
    i = 1;
    R = 188.92; % Ideal gas Constant
    temp(1) = t_init;
    rho1 = refpropm('D','T',t_init,'P',current_pressure/1e3, 'CO2');
    
    u = refpropm('U','T',t_init,'P',current_pressure/1e3, 'CO2');
    
    hgl1 = current_pressure(1)/(rho1 * gravity) + 3200;
    direction = 1;
    
    pressure_loss = 0;
    delta_l = 1;
    gravity_gain = 0;
    % Start at top
    % when each range covers h  --> h + delta_l
    % 0 definited as top surface
    % 3200 defined as underground
    i = 2;

   
    
    for height = 0 : delta_l : 3200 - delta_l
        
        temp(i) = temp(i-1);
        init_guess = temp(i);
        output_temp = 0;
        while abs(init_guess - output_temp > .01)
            init_guess = temp(i);
            dynamic_v1 = refpropm('V','T',temp(i),'P',current_pressure(i-1)/1e3, 'CO2');
            rho1(i) = refpropm('D','T',temp(i),'P',current_pressure(i-1)/1e3, 'CO2');
            velo1 = m_dot/(rho1(i) * area);
            reynolds(i) = rho1(i) * velo1 * diameter/dynamic_v1;

            if reynolds(i) >= 4000
                f1 = -1.8 * log10((6.9 / reynolds(i)) + (e_d_ratio/ 3.7)^1.); % turb only
                f(i) = (1/f1)^2;
            else
                f(i) = 64/reynolds(i);
            end
            
            head_loss(i) = (f(i) * delta_l * velo1^2) / (diameter * 2  * gravity);
            pressure_drop(i) = rho1(i) * gravity * head_loss(i); %pa
            
            u(i) = gravity * head_loss(i) + u(i-1);
            
            current_pressure(i) = current_pressure(i-1) - pressure_drop(i);
            rho2 = refpropm('D','T',temp(i),'P',current_pressure(i)/1e3, 'CO2');
            velo2 = m_dot/(rho2 * area);
            velo_diff = (velo2^2 - velo1^2)/2;
            
            % Gravity, negative delta l since 2 is below 1
            
            % Gravity, negative delta l since 2 is below 1
            % Change to P/Rho Z * G
            % Average p and rho to get Z
            % Z in refprop
            Z = refpropm('Z','T',temp(i),'P',current_pressure(i)/1e3, 'CO2');
            comp_const = current_pressure(i-1)/(((rho2 + rho1(i))/2)* Z);            
            %p2 = current_pressure(i-1) * exp((-1/(R * temp(i))) * (velo_diff +  gravity * -1 * delta_l));
            p2 = current_pressure(i-1) * exp((-1/(comp_const)) * (velo_diff +  gravity * -1 * delta_l));
            
            %p2 = current_pressure(i-1) * exp((-1/(R * temp(i))) * (velo_diff +  gravity * -1 * delta_l));
            gg = p2 - current_pressure(i-1);
            
            %cp2(i) = cp2(i-1) +  rho1(i) * ((velo1^2 - velo2^2) * 1/2 + gravity * delta_l) - pressure_drop(i); 
            current_pressure(i) = current_pressure(i-1) + gg -  pressure_drop(i);
            
            temp(i-1) = refpropm('T','P',current_pressure(i)/1e3,'U',u(i), 'CO2');
            
            output_temp = HT(m_dot, current_pressure(i-1),temp(i-1), delta_l, direction, current_pressure(i), height);
            temp(i) = output_temp;
            u(i) = refpropm('U','T',temp(i),'P',current_pressure(i)/1e3, 'CO2');
        end
        
        gravity_gain = gravity_gain + gg; 
        pressure_loss = pressure_loss + pressure_drop(i);
           
        i = i + 1;
    end
    current_press = current_pressure(end);
    
end

