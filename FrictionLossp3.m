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

temp = (15 + 129)/2 + 273;
area = pi * (diameter/2)^2 ;
md = 3.24:.1:18;
i = 1

leg =20:10:80%pascals, can set this one
% p3_val = 10e6:20e6:100e6
i = 1
k = 1e-12 ;% 
A = 20 ;% m^-1
L = 112 ;% m
Pa = 32.5e6 ;
Tb = 129 + 273; % Temp at bottom
p3 = 20e6 :10e6:80e6;
md = 3.24 : 1 : 15.24;
%md = 10
%for p3  = 20e6 :10e6:80e6
md = 5
for p3  = 40e6
    for j = 1:length(md)
        [p4(i,j), friction_loss_down(i,j), gravity_gain_down(i,j), temp_down] = pressure_drop_down(md(j), p3);       
        density(i,j) = refpropm('D','T',Tb,'P',p4(i,j)/ 1e3, 'CO2');
        dynamic_v(i,j) = refpropm('V','T',Tb,'P',p4(i,j)/1e3, 'CO2');
      
        final_mdot(j) = md(j) - 3.24;
        Q(i,j) = final_mdot(j) / density(i,j);
        delta_p(i,j) = Q(i,j) * dynamic_v(i,j) * L/(k*A);
        
        p5(i,j) = p4(i,j) - delta_p(i,j);
        [p6(i,j), friction_loss_up(i,j), gravity_loss_up(i,j), temp_top(i,j), temp_up, gg, current_p] = pressure_drop_up(final_mdot(j), p5(i,j));   
        final_delta_p(i,j) = p3 - p6(i,j) 
        
    end 
    i = i+1;
end
%%
for iN = 1:length(leg)
      legendCell{iN} = num2str(leg(iN),'P3=%-d');
end
 
%%
% for ind = 2:length(temp_down)
%     dt(ind-1) = temp_down(ind) - temp_down(ind-1);
%     dt2(ind -1) = temp_up(ind) - temp_up(ind-1);
% end
% delta_l = 1;
% figure
% plot(0 : delta_l : 3200 - delta_l, dt)
% yyaxis right
% plot(0 : delta_l : 3200, temp_down)
% title('t down')
% 
% 
% figure
% plot(0 : delta_l : 3200 - delta_l, dt2)
% yyaxis right
% plot(0 : delta_l : 3200, temp_up)







%%
% figure
% plot(md, head_loss_up)
% title('Friction head loss due to friction in a pipe, up ')
% xlabel('Mass Flow Rate (kg/s)')
% ylabel('Friction Head (m)')
% legend(legendCell)


%legend('10' ,'20', '30', '40' ,'50')

% figure
% hold on 
% plot(final_mdot - 3.24, gravity_loss_up)
% title('Gravity loss up the pipe')
% xlabel('Mass Flow Rate Down(kg/s)')
% ylabel('Pressure loss (Pa)')
% %legend('10' ,'20', '30', '40' ,'50')
% legend(legendCell)
% %%
% figure
% hold on 
% plot(md, gravity_gain_down)
% title('Gravity gain down the pipe')
% xlabel('Mass Flow Rate Down(kg/s)')
% ylabel('Pressure gain (Pa)')
% %legend('10' ,'20', '30', '40' ,'50')
% legend(legendCell)
% %%
% figure
% hold on 
% plot(md, gravity_gain_down - gravity_loss_up)
% title('Gravity gain (difference)')
% xlabel('Mass Flow Rate Down(kg/s)')
% ylabel('Pressure gain (Pa)')
% %legend('10' ,'20', '30', '40' ,'50')
% legend(legendCell)
% 
% 
% 
% %%
% figure
% hold on 
% plot(md, press_drop_down)
% title('Pressure loss due to friction in a pipe, down')
% xlabel('Mass Flow Rate Down(kg/s)')
% ylabel('Pressure Drop (Pa)')
% %legend('10' ,'20', '30', '40' ,'50')
% legend(legendCell)
% 
% 
% figure
% hold on 
% plot(md, press_drop_up + press_drop_down)
% title('Pressure loss due to friction in a pipe, both')
% xlabel('Mass Flow Rate Down(kg/s)')
% ylabel('Pressure Drop (Pa)')
% legend(legendCell)
% 
% figure
% hold on 
% plot(md, press_drop_down)
% title('Pressure loss due to friction in a pipe, up')
% xlabel('Mass Flow Rate Down(kg/s)')
% ylabel('Pressure Drop (Pa)')
% legend(legendCell)
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
% figure
% plot(md,p6) 
% xlabel('Mass Flow Rate, Down')
% ylabel('Pressure 6 (pa)')
% title('P6 vs Mass Flow Rate for various P3')
% legend(legendCell)
% 
% %%
figure
plot(md, final_delta_p) 
xlabel('Mass Flow Rate, Down')
ylabel('Pressure Difference (3-6) (pa)')
title('Final Delta P vs Mass Flow Rate for various P3')
legend(legendCell)
% 
% %%
% figure
% plot(md, p4)
% xlabel('Mass Flow Rate, Down')
% ylabel('Pressure 4 (pa)')
% title('P4 vs Mass Flow Rate for various P3')
% legend(legendCell)
% 
% figure
% plot(md, delta_p)
% xlabel('Mass Flow Rate, Down')
% ylabel('Delta P (pa)')
% title('Delta P Accross Sand vs Mass Flow Rate for various P3')
% legend(legendCell)
% 
% figure
% plot(md, p5)
% xlabel('Mass Flow Rate, Down')
% ylabel('P5 (pa)')
% title('P5 vs Mass Flow Rate for various P3')
% legend(legendCell)
% %%
% figure
% plot(md, density)
% xlabel('Mass Flow Rate, Down')
% ylabel('Density')
% title('Density vs Mass Flow Rate for various P3')
% legend(legendCell)
% %%
% figure
% plot(md, dynamic_v)
% xlabel('Mass Flow Rate, Down')
% ylabel('Dyn V')
% title('Dyn V vs Mass Flow Rate for various P3')
% legend(legendCell)
% 
% figure
% plot(md, Q)
% xlabel('Mass Flow Rate, Down')
% ylabel('Vol Flow Rate')
% title('Vol Flow Rate vs Mass Flow Rate for various P3')
% legend(legendCell)

