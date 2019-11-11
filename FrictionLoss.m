m_dot = 3.4; 
diameter = 2.4; % inch
diameter = 0.06096; % meter
pressure = 6; % Mpa
pressure = 6e6
temp = 15 + 273;
rho = refpropm('D','T',temp,'P', pressure /1e3, 'CO2');
area = pi * (diameter/2)^2 ;
velo = m_dot /(rho * area);
dynamic_v = refpropm('V','T',temp,'P',pressure/1e3, 'CO2');
reynolds = rho * velo * diameter/dynamic_v ;

roughness = .05 % mm, from https://neutrium.net/fluid_flow/absolute-roughness/
% new = .02 -.05 mm, slightly corroded = .05 - .15, moderate corrosion =
% .15-1mm

e_d_ratio = roughness/1000 * diameter

f1 = -1.8 * log10((6.9 / reynolds) + (e_d_ratio/ 3.7)^1.)
f = (1/f1)^2

length = 3.2 * 1000 % meters
gravity = 9.81

head_loss = (f * length * velo^2) / (diameter * 2  * gravity)
pressure_drop = rho * gravity * head_loss %pa
