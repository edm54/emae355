close all
clear all

set(0, 'DefaultAxesFontWeight', 'normal', ...
      'DefaultAxesFontSize', 18, ...
      'DefaultAxesFontAngle', 'normal', ... 
      'DefaultAxesFontWeight', 'normal', ... 
      'DefaultAxesTitleFontWeight', 'bold', ...
      'DefaultAxesTitleFontSizeMultiplier', 1.2) ;
set(groot,'defaultLineLineWidth',3)


roughness = .05 ;% mm, from https://neutrium.net/fluid_flow/absolute-roughness/
% new = .02 -.05 mm, slightly corroded = .05 - .15, moderate corrosion =
% .15-1mm

L1 = 3.2 * 1000; % meters
gravity = 9.81;

diameter = 2.4; % inch
diameter = 0.06096; % meter
e_d_ratio = roughness/1000 * diameter;
pressure = 10e6 %Pa
temp = (15 + 129)/2 + 273;
rho = refpropm('D','T',temp,'P', pressure /1e3, 'CO2');
area = pi * (diameter/2)^2 ;
md = 3.24:.1:12
i = 1
% up 
% for m_dot = 3.24:.1:17
%     velo_up = (m_dot - 3.24) /(rho * area);
%     velo_down = (m_dot) /(rho * area);
%     
%     dynamic_v = refpropm('V','T',temp,'P',pressure/1e3, 'CO2');
%     
%     reynolds_up = rho * velo_up * diameter/dynamic_v ;
%     reynolds_down = rho * velo_down * diameter/dynamic_v ;
% 
%     f1_up = -1.8 * log10((6.9 / reynolds_up) + (e_d_ratio/ 3.7)^1.);
%     f_up(i) = (1/f1_up)^2;
% 
%     f1_down = -1.8 * log10((6.9 / reynolds_down) + (e_d_ratio/ 3.7)^1.);
%     f_down(i) = (1/f1_down)^2;
%     
%     head_loss_up(i) = (f_up(i) * L1 * velo_up^2) / (diameter * 2  * gravity);
%     pressure_drop_up(i) = rho * gravity * head_loss_up(i); %pa
%     
%     head_loss_down(i) = (f_down(i) * L1 * velo_down^2) / (diameter * 2  * gravity);
%     pressure_drop_down(i) = rho * gravity * head_loss_down(i); %pa
%     
%     i = i+1
% 
% end
% %%
% 
% figure
% plot(md, head_loss_up)
% title('Friction head loss due to friction in a pipe, up ')
% xlabel('Mass Flow Rate (kg/s)')
% ylabel('Friction Head (m)')
% 
% figure
% hold on 
% plot(md, pressure_drop_up)
% plot(md, pressure_drop_down)
% plot(md, pressure_drop_up + pressure_drop_down)
% title('Pressure loss due to friction in a pipe')
% xlabel('Mass Flow Rate Down(kg/s)')
% ylabel('Pressure Drop (Pa)')
% legend('Up', 'Down', 'Both')
% 
% figure
% plot(md , head_loss_down)
% title('Friction head loss due to friction in a pipe, down ')
% xlabel('Mass Flow Rate (kg/s)')
% ylabel('Friction Head (m)')
%% 
leg =50:5:70%pascals, can set this one
p6_val = 50e6:5e6:70e6
i = 1
k = 1e-12 ;% 
A = 20 ;% m^-1
L = 112 ;% m
Pa = 32.5e6 ;
Tb = 129 + 273; % Temp at bottom
for p6  = 50e6:5e6:70e6
    for j = 1:length(md)
%         dynamic_v_top = refpropm('V','T',temp,'P',p6/1e3, 'CO2');
%         rho(i,j) = refpropm('D','T',temp,'P',p6/1e3, 'CO2');
%         %rho(i,j) = 400
%         velo_up = (md(j) - 3.24) /(rho(i,j) * area);
%         velo_down = md(j)/(rho(i,j) * area);
%         reynolds_up = rho(i,j) * velo_up * diameter/dynamic_v_top ;
%         reynolds_down = rho(i,j) * velo_down * diameter/dynamic_v_top ;
% 
%         f1_up = -1.8 * log10((6.9 / reynolds_up) + (e_d_ratio/ 3.7)^1.);
%         f_up(i,j) = (1/f1_up)^2;
% 
%         f1_down = -1.8 * log10((6.9 / reynolds_down) + (e_d_ratio/ 3.7)^1.);
%         f_down(i,j) = (1/f1_down)^2;
% 
%         head_loss_up(i,j) = (f_up(i,j) * L1 * velo_up^2) / (diameter * 2  * gravity);
%         pressure_drop_up(i,j) = rho(i,j) * gravity * head_loss_up(i,j); %pa
% 
%         head_loss_down(i,j) = (f_down(i,j) * L1 * velo_down^2) / (diameter * 2  * gravity);
%         pressure_drop_down(i,j) = rho(i,j) * gravity * head_loss_down(i,j); %pa      
%         
        
        press_drop_up(i,j) = pressure_drop_up(md(j)-3.24, p6);
        p5(i,j) = p6 + press_drop_up(i,j);
        
        density(i,j) = refpropm('D','T',Tb,'P',p5(i,j)/ 1e3, 'CO2');
        dynamic_v(i,j) = refpropm('V','T',Tb,'P',p5(i,j)/1e3, 'CO2');
      
        final_mdot(j) = md(j) - 3.24;
        Q(i,j) = final_mdot(j) / density(i,j);
        delta_p(i,j) = Q(i,j) * dynamic_v(i,j) * L/(k*A);
        % + pa
        p4(i,j) = p5(i,j) + delta_p(i,j);
        gravity_gain = 0; % TODO
        press_drop_down(i,j) = pressure_drop(md(j), p4(i,j));
        p3(i,j) = p4(i,j) + press_drop_down(i,j) + gravity_gain;
        delta_p(i,j) = p3(i,j) - p6;
        
    end
%     for j = 1:length(md)
%       
%         p5(i,j) = p6 + pressure_drop_up(i,j);
%         
%         density(i,j) = refpropm('D','T',Tb,'P',p5(i,j)/ 1e3, 'CO2');
%         dynamic_v(i,j) = refpropm('V','T',Tb,'P',p5(i,j)/1e3, 'CO2');
%       
%         final_mdot(j) = md(j) - 3.24;
%         Q(i,j) = final_mdot(j) / density(i,j);
%         delta_p(i,j) = Q(i,j) * dynamic_v(i,j) * L/(k*A);
%         p4(i,j) = Pa + delta_p(i,j);
%         gravity_gain = 0; % TODO
%         p3(i,j) = p4(i,j) + pressure_drop_down(i,j) + gravity_gain;
%     end
    i = i+1
end
%%
for iN = 1:length(leg)
      legendCell{iN} = num2str(leg(iN),'P6=%-d');
end
 
%%
% figure
% plot(md, head_loss_up)
% title('Friction head loss due to friction in a pipe, up ')
% xlabel('Mass Flow Rate (kg/s)')
% ylabel('Friction Head (m)')
% legend(legendCell)

figure
hold on 
plot(md, press_drop_up)
title('Pressure loss due to friction in a pipe, up')
xlabel('Mass Flow Rate Down(kg/s)')
ylabel('Pressure Drop (Pa)')
legend(legendCell)
%legend('10' ,'20', '30', '40' ,'50')




figure
hold on 
plot(md, press_drop_down)
title('Pressure loss due to friction in a pipe, down')
xlabel('Mass Flow Rate Down(kg/s)')
ylabel('Pressure Drop (Pa)')
%legend('10' ,'20', '30', '40' ,'50')
legend(legendCell)

%%
figure
hold on
plot(md, press_drop_down(2, :))
legend('old', 'new')
% plot(md, pressure_drop_down(2, :))
% plot(md, pressure_drop_down1(2, :))
%%
%legend(legendCell)

figure
hold on 
plot(md, press_drop_up + press_drop_down)
title('Pressure loss due to friction in a pipe, both')
xlabel('Mass Flow Rate Down(kg/s)')
ylabel('Pressure Drop (Pa)')
legend(legendCell)

% figure
% hold on 
% plot(md, pressure_drop_up)
% plot(md, pressure_drop_down)
% plot(md, pressure_drop_up + pressure_drop_down)
% title('Pressure loss due to friction in a pipe')
% xlabel('Mass Flow Rate Down(kg/s)')
% ylabel('Pressure Drop (Pa)')
% legend('Up', 'Down', 'Both')



% figure
% plot(md , head_loss_down)
% title('Friction head loss due to friction in a pipe, down ')
% xlabel('Mass Flow Rate (kg/s)')
% ylabel('Friction Head (m)')
% legend(legendCell)
%%
figure
plot(md , rho)
title('Rho')
xlabel('Mass Flow Rate (kg/s)')
ylabel('rho kg')
legend(legendCell)



%%
figure
plot(md,p3) 
xlabel('Mass Flow Rate, Down')
ylabel('Pressure 3 (pa)')
title('P3 vs Mass Flow Rate for various P6')
legend(legendCell)

%%
figure
plot(md, delta_p) 
xlabel('Mass Flow Rate, Down')
ylabel('Pressure Difference (3-6) (pa)')
title('Delta P vs Mass Flow Rate for various P6')
legend(legendCell)

%%
figure
plot(md, p4)
xlabel('Mass Flow Rate, Down')
ylabel('Pressure 4 (pa)')
title('P4 vs Mass Flow Rate for various P6')
legend(legendCell)

figure
plot(md, delta_p)
xlabel('Mass Flow Rate, Down')
ylabel('Pelta P (pa)')
title('Delta P vs Mass Flow Rate for various P6')
legend(legendCell)

figure
plot(md, p5)
xlabel('Mass Flow Rate, Down')
ylabel('P5 (pa)')
title('P5 vs Mass Flow Rate for various P6')
legend(legendCell)
%%
figure
plot(md, density)
xlabel('Mass Flow Rate, Down')
ylabel('Density')
title('Density vs Mass Flow Rate for various P6')
legend(legendCell)
%%
figure
plot(md, dynamic_v)
xlabel('Mass Flow Rate, Down')
ylabel('Dyn V')
title('Dyn V vs Mass Flow Rate for various P6')
legend(legendCell)

figure
plot(md, Q)
xlabel('Mass Flow Rate, Down')
ylabel('Vol Flow Rate')
title('Vol Flow Rate vs Mass Flow Rate for various P6')
legend(legendCell)

