function [cp, pressure_loss, gravity_gain_total, temp_final] = pressure_drop_up(m_dot, pressure1)
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
    
    i = 1;
    R = 188.92; % Ideal gas Constant
    %rho1 = refpropm('D','T',temp,'P',current_pressure/1e3, 'CO2');
    pressure_loss = 0;
    delta_l = 1;
    t_bottom = 129 + 273;
    temp(1) = t_bottom;
    i = 2;
    direction = 2;
    gravity_gain = 0;
    % Start at bottom
    % when each range covers h  --> h + delta_l
    temp_final = 0;
    
    if m_dot > 0
        for height = 3200 - delta_l : -1 * delta_l : 0 

            temp(i) = temp(i-1);
            init_guess = temp(i);
            output_temp = 0;
            while abs(init_guess - output_temp) > .01
                init_guess = temp(i);
                dynamic_v1 = refpropm('V','T',temp(i),'P',current_pressure(i-1)/1e3, 'CO2');
                rho1(i) = refpropm('D','T',temp(i),'P',current_pressure(i-1)/1e3, 'CO2');

                velo1 = m_dot/(rho1(i) * area);
                velo(i) = m_dot/(rho1(i) * area);
                reynolds(i) = rho1(i) * velo1 * diameter/dynamic_v1;

                if reynolds(i) >= 4000
                    f1 = -1.8 * log10((6.9 / reynolds(i)) + (e_d_ratio/ 3.7)^1.); % turb only
                    f(i) = (1/f1)^2;
                else
                    f(i) = 64/reynolds(i);
                end

                head_loss(i) = (f(i) * delta_l * velo1^2) / (diameter * 2  * gravity);
                pressure_drop(i) = rho1(i) * gravity * head_loss(i); %pa
                
                current_pressure(i) = current_pressure(i-1) - pressure_drop(i);
                rho2 = refpropm('D','T',temp(i),'P',current_pressure(i)/1e3, 'CO2');
                velo2 = m_dot/(rho2 * area);

                velo_diff = (velo2^2 - velo1^2)/2;
                % Gravity
                
                Z = refpropm('Z','T',temp(i),'P',current_pressure(i)/1e3, 'CO2');
                comp_const = current_pressure(i-1)/(((rho2 + rho1(i))/2)* Z) ;           
                %p2 = current_pressure(i-1) * exp((-1/(R * temp(i))) * (velo_diff +  gravity *  delta_l));
                p2 = current_pressure(i-1) * exp((-1/(comp_const)) * (velo_diff +  gravity * delta_l));
            
                 %p2 = current_pressure(i-1) * exp((-1/(R * temp(i))) * (velo_diff +  gravity * delta_l));
                 gg(i) = current_pressure(i-1) - p2; 
                 gravity_gain(i) = gravity_gain(i-1) + current_pressure(i-1) - p2;
                 %current_pressure(i) = p2;
                 current_pressure(i) = current_pressure(i-1) - gg(i) - pressure_drop(i);
                 
                 output_temp = HT(m_dot, current_pressure(i-1),temp(i-1), delta_l, direction, current_pressure(i), height);
                 temp(i) = output_temp;
                 
            end
            pressure_loss = pressure_loss + pressure_drop(i);
            

             i = i + 1;
        end
    end
    cp = current_pressure(end);
    gravity_gain_total = gravity_gain(end);
    temp_final = temp(end);
end

