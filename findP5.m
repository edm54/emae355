function [p5, t5, p3, t3, p4, t4] = findP5(p2,md)

k = 1e-12 ;% 
A = 20 ;% m^-1
L = 112 ;% m
Pa = 32.5e6 ;
Tb = 129 + 273; % Temp at bottom

[p3, friction_loss_down, gravity_gain_down, temp_down] = pressure_drop_down(md, p2);       
density = refpropm('D','T',Tb,'P',p3/ 1e3, 'CO2');
dynamic_v = refpropm('V','T',Tb,'P',p3/1e3, 'CO2');

t3 = temp_down;

final_mdot = md - 3.24;
Q = final_mdot/ density;
delta_p = Q * dynamic_v* L/(k*A);
t4 = 129+273;
p4 = p3 - delta_p;
[p5, friction_loss_up, gravity_loss_up, temp_top] = pressure_drop_up(final_mdot, p4);   
t5 = temp_top;
        
end

