close all
clear all

roughness = .05 ;% mm, from https://neutrium.net/fluid_flow/absolute-roughness/
% new = .02 -.05 mm, slightly corroded = .05 - .15, moderate corrosion =
% .15-1mm

e_d_ratio = roughness/1000 * diameter;

diameter = 2.4; % inch
diameter = 0.06096; % meter
pressure = 6e6 %Pa
temp = 15 + 273;
rho = refpropm('D','T',temp,'P', pressure /1e3, 'CO2');
area = pi * (diameter/2)^2 ;
md = 3.24:.1:13.24
i = 1
% up 
for m_dot = 3.24:.1:13.24
    velo_up = (m_dot - 3.24) /(rho * area);
    dynamic_v = refpropm('V','T',temp,'P',pressure/1e3, 'CO2');
    reynolds_up = rho * velo_up * diameter/dynamic_v ;


    f1 = -1.8 * log10((6.9 / reynolds) + (e_d_ratio/ 3.7)^1.);
    f(i) = (1/f1)^2;

    length = 3.2 * 1000; % meters
    gravity = 9.81;

    head_loss(i) = (f(i) * length * velo^2) / (diameter * 2  * gravity);
    pressure_drop1(i) = rho * gravity * head_loss(i); %pa
    i = i+1

end
%%


figure
plot(md , head_loss)
title('Friction head loss due to friction in a pipe, up ')
xlabel('Mass Flow Rate (kg/s)')
ylabel('Friction Head (m)')
figure
plot(md, pressure_drop1)
hold on 

md = 3.24:.1:13.24
i = 1
% Down the pipe
for m_dot = 3.24:.1:13.24
    velo = (m_dot) /(rho * area);
    dynamic_v = refpropm('V','T',temp,'P',pressure/1e3, 'CO2');
    reynolds = rho * velo * diameter/dynamic_v ;

    roughness = .05 ;% mm, from https://neutrium.net/fluid_flow/absolute-roughness/
    % new = .02 -.05 mm, slightly corroded = .05 - .15, moderate corrosion =
    % .15-1mm

    e_d_ratio = roughness/1000 * diameter;

    f1 = -1.8 * log10((6.9 / reynolds) + (e_d_ratio/ 3.7)^1.);
    f(i) = (1/f1)^2;

    length = 3.2 * 1000; % meters
    gravity = 9.81;

    head_loss(i) = (f(i) * length * velo^2) / (diameter * 2  * gravity);
    pressure_drop2(i) = rho * gravity * head_loss(i); %pa
    i = i+1

end

%%
figure
hold on 
plot(md, pressure_drop1)
plot(md, pressure_drop2)
plot(md, pressure_drop1 + pressure_drop2)
title('Pressure loss due to friction in a pipe')
xlabel('Mass Flow Rate Down(kg/s)')
ylabel('Pressure Drop (Pa)')
legend('Up', 'Down', 'Both')

figure
plot(md , head_loss)
title('Friction head loss due to friction in a pipe, down ')
xlabel('Mass Flow Rate (kg/s)')
ylabel('Friction Head (m)')
%% 

p6 = 12e6 %pascals, can set this one
p5 = p6 - pressure_drop
