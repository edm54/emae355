clear all
clc
close all
%[q, t] = radiation(306, 295, 420, 1.5)

global alp kin_v  Pr k g B
g = 9.81; %m/s 

Ts = 160 + 273;
Ti = 25 + 273 ;%80 f, room temperature
T_average = (Ts + Ti)/2;
B = 1/T_average;
% Prandtl number
Pr = refpropm('^','T',Ti,'P',101, 'Air.ppf');
% V   Dynamic viscosity [Pa*s], Mu
Dynamic_V = refpropm('V','T',Ti,'P',101, 'Air.ppf');

%  $   Kinematic viscosity [cm^2/s], Nu
kin_v = refpropm('$','T',Ti,'P',101, 'Air.ppf');
kin_v = kin_v /100^2 ;
kin_v = 26.4e-6;
%Kin_v = Dynamic_v/rho

%   Thermal diffusivity [cm^2/s], alpha
alp = refpropm('%','T',Ti,'P',101, 'Air.ppf');
alp = alp/ 100^2;
alp = 38.3e-6;

%L   Thermal conductivity [W/(m K)]
k = refpropm('L','T',Ti,'P',101, 'Air.ppf')

faces = 4 ;
t_total = 7 * 60; % seconds
%k = .525 % w/m from us department of commerce, cast iron 

%k = 80 % w/mk from a table carter found, young ch 15 

%% Inout Paramaters
%exp((velo2^2 - velo1^2)/(2*R*(T1)))

area_flat= [.891]; % M2 for ONE FLAT Face
perimeter_flat = 3.82; %M
l_flat = area_flat./ [perimeter_flat];

% Total area of 4 vertical walls
area_vert = [.65/2]; %meter2;
perimeter_vert = 8.32/2;
l_vert = [perimeter_vert]./ area_vert;

%% Vertical Walls
t_vert = 7*60; %sec
[q_vertical, Q_vertical] = vertical_wall_convection(Ts, Ti, t_vert, l_vert,  area_vert)

%% Convection from bottom
t_bottom_convection= 3*60;
[q_bottom, Q_bottom] = free_convection_up(Ts, Ti, t_bottom_convection, l_flat, area_flat)

%% Floating
t_float_t= 30; % seconds

[q_float_b, Q_float_b] =  free_convection_down(Ts, Ti, t_float_t, l_flat, area_flat)
[q_float_t, Q_float_t] = free_convection_up(Ts, Ti, t_float_t, l_flat, area_flat)

q_float = q_float_b  + q_float_t;
Q_float = Q_float_b + Q_float_t

%% top plate convection
t_open_t= 2.5*60; % seconds
[q_top, Q_top] = free_convection_down(Ts, Ti, t_open_t, l_flat, area_flat);

%% Forced Convection
t_forced_t= 45; % naturally convect same 
[q_forced, Q_forced] = forced_convection(Ts, Ti, t_forced_t, .1*l_flat, .1*area_flat)

%% Radiation
t_rad = 7 * 60;
[q_rad, Q_rad] = radiation(Ts, Ti, t_rad, 2.*area_vert);

%% New Calcs
% 
t_rad_n = 6 * 60
t_vert_n = 6 * 60
t_bottom = 2 * 60
t_top = 1.5 * 60

[q_rad_n, Q_rad_n] = radiation(Ts, Ti, t_rad_n , 2.*area_vert);
[q_vertical_n, Q_vertical_n] = vertical_wall_convection(Ts, Ti, t_vert_n,  l_vert, area_vert)
[q_top_n, Q_top_n] = free_convection_down(Ts, Ti, t_top, l_flat, area_flat);
[q_bottom_n, Q_bottom_n] = free_convection_up(Ts, Ti, t_bottom, l_flat, area_flat)


%%
Q_absolute_total = Q_forced + Q_bottom + Q_float ...
                 + Q_top + 2*Q_vertical + Q_rad

Q_new_total = Q_forced  + Q_bottom_n + Q_float ...
                 + Q_top_n + 2*Q_vertical_n + Q_rad_n

savings = (Q_absolute_total - Q_new_total)/Q_absolute_total
     
Q_loss = Q_absolute_total - Q_new_total

figure
pie([Q_forced Q_bottom Q_float ...
                 Q_top  2*Q_vertical  Q_rad])
legend('Forced', 'Bottom' , 'Float', 'Top', 'Vertical', 'Radiation')
title('Original Heat Transfer')

figure
pie([Q_forced Q_bottom_n Q_float ...
                 Q_top_n  2*Q_vertical_n  Q_rad_n Q_loss])
legend('Forced', 'Bottom' , 'Float', 'Top', 'Vertical', 'Radiation', 'Savings')
title('New Heat Transfer')


figure
pie([q_forced q_bottom q_float ...
                 q_top  2*q_vertical  q_rad])
legend('Forced', 'Bottom' , 'Float', 'Top', 'Vertical', 'Radiation')
title('q')

%%
total_toys_per_machine = 24 * 60 /7
time_with_new_cycle = total_toys_per_machine * 6/60
time_saved = 24 - time_with_new_cycle
% Min number of machines to turn off one
min_machines = 24/time_saved + 1

% Heatloss per cycle * cylces per machine * machines
heat_loss_daily_original = Q_absolute_total * total_toys_per_machine * 8
heat_loss_daily_new = Q_new_total * total_toys_per_machine * 7

total_savings = (heat_loss_daily_original - heat_loss_daily_new)/heat_loss_daily_original


function[q_vert, Q_total ] = vertical_wall_convection(Ts,Ti,time, l_vert, area_vert)
global alp kin_v Pr k g B
% Vertical walls
Ra_vert = (g * B * (Ts - Ti).*l_vert.^3)/(alp * kin_v);
nusselt_vert = (.825 + (.387.*Ra_vert.^(1/6))/(1 + (.492/Pr)^(9/16))^(8/27)).^2;

h_vert = nusselt_vert * k./l_vert
q_vert = h_vert .* area_vert * (Ts - Ti); % watts 
% copy code, for diff length
Q_total = q_vert * time ;

end

function [q_forced, Q_total] = forced_convection(Ts, Ti, time, L_forced, area_forced)
velo2=200; %m/s
global kin_v Pr k 
%L_forced = .05 %char len
%area_forced = .025

% Recalculate kin_V, Pr
Rex = velo2 * L_forced/kin_v;

nusselt_forced = .0296 * Rex ^.8 * Pr ^(1/3);
h_forced = nusselt_forced * k./L_forced;
q_forced = h_forced .* area_forced * (Ts - Ti);
Q_total = q_forced * time;
end


function [q_down, Q_total] = free_convection_down(Ts, Ti, time, L, area_flat)
%Facing down
%same L and A as bottom  plate
global alp kin_v k g B

Ra_down = (g * B * (Ts - Ti).*L.^3)/(alp * kin_v)
nusselt_down = .52 * Ra_down^(1/5)
h_down = nusselt_down * k./L
q_down = h_down .* area_flat * (Ts - Ti); % watts 
Q_total = q_down * time;
end


function [q_rad, Q_total] = radiation(Ts, Ti, time, area_vert)
emiss= .3;
boltzman= 5.67e-8; %W/(m^2K^4)
q_rad= emiss*boltzman*(Ts^4-Ti^4)*area_vert;
Q_total  = q_rad * time;
end

function [q_free_up, Q_total] = free_convection_up(Ts, Ti, time, l_flat, area_flat)
% Flt plate facing up
global g B alp kin_v k
Ra_flat = (g * B * (Ts - Ti).*l_flat.^3)/(alp * kin_v);
nusselt_flat = .15 * Ra_flat.^(1/3);
h_up = nusselt_flat * k./l_flat
q_free_up = h_up .* area_flat * (Ts - Ti); % watts 
Q_total = q_free_up * time;

end



%% Nozzle

% T1 = 20+273 %C
% P1 = 620.528 %kps
% D1 = .00953 % m
% D2 = 0.00635 %m  from 1/4 inc
% 
% rho1 = refpropm('D','T',T1,'P',P1, 'Air.ppf');
% 
% v_rate1 = .0165 %m3/sec
% A1 = pi * (D1/2)^2
% A2 = pi * (D2/2)^2
% velo1 = v_rate1/A1
% 
% velo2 = min(velo1 * A1/A2, 340)
% velo2 = 340
% gamma = 1.4
% Tstag = T1*(gamma + 1)/2
% Pstag = P1 * ((gamma + 1)/2)^(gamma/(gamma - 1))
% P2 = Pstag/(((gamma - 1)/gamma)+1)
% 
% T2 = Tstag / (((gamma-1)/2) + 1)
% rho2 = (rho1 * velo1 * A1)/ (velo2*A2)
% 
% 
% R = 8.314 %j/molk
% 
% P2 = P1*exp((velo1^2 - velo2^2)/(2*R*(T1+273)))
% 
% %P2 = P1 + rho1 * (velo1^2)/2  -. rho2 * (velo2^2)/2


