global alp kin_v  Pr k g B
g = 9.81 %m/s 
T_average = (Ts + Ti)/2
B = 1/T_average
Ts = 160 + 273
Ti = 26.6 + 273 %80 f, room temperature

% Prandtl number
Pr = refpropm('^','T',Ti,'P',101, 'Air.ppf');
% V   Dynamic viscosity [Pa*s], Mu
Dynamic_V = refpropm('V','T',Ti,'P',101, 'Air.ppf');

%  $   Kinematic viscosity [cm^2/s], Nu
kin_v = refpropm('$','T',Ti,'P',101, 'Air.ppf');
kin_v = kin_v /100^2 
kin_v = 26.4e-6
%Kin_v = Dynamic_v/rho

%   Thermal diffusivity [cm^2/s], alpha
alp = refpropm('%','T',Ti,'P',101, 'Air.ppf');
alp = alp/ 100^2
alp = 38.3e-6

%L   Thermal conductivity [W/(m K)]
k = refpropm('L','T',Ti,'P',101, 'Air.ppf');

faces = 4 
t_total = 7 * 60 % seconds
%k = .525 % w/m from us department of commerce, cast iron 

%k = 80 % w/mk from a table carter found, young ch 15 



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


%% Total Heat Transfer over time
%exp((velo2^2 - velo1^2)/(2*R*(T1)))
l_flat = [.75] % todo: L should be As/perimeter
area_flat= [.625]

l_vert = [.5 .25]
area_vert = [ .25, .125] %meter2

%% Closed
t_closed= 4*60; %sec
q_closed = 0
Q_total_closed = 0
q_closed, Q_total_closed = vertical_wall_convection(Ts, Ti, t_closed, l_vert, area_vert)

%% Open
t_open_b= 3*60;
q_bottom = 0
Q_total_bottom = 0
q_bottom, Q_total_bottom = free_convection_down(Ts, Ti, t_open_b, l_flat, area_flat)

%% six sides are open
t_float_t= 30; 

q_bottom_all = 0
Q_bottom_all = 0
q_closed_all = 0
Q_closed_all = 0
q_top_all = 0
Q_top_all = 0

q_bottom_all, Q_bottom_all =  free_convection_down(Ts, Ti, t_float_t, l_flat, area_flat)
q_closed_all, Q_closed_all = vertical_wall_convection(Ts, Ti, t_float_t, l_vert, area_vert)
q_top_all, Q_top_all = free_convection_up(Ts, Ti, t_float_t, l_flat, area_flat)

q_all_sides = q_bottom_all + q_closed_all + q_top_all
Q_total_all_sides = Q_bottom_all + Q_closed_all + Q_top_all

%% 1 face and 4 vertical walls
t_open_t= 2.5*60; % 
q_vert = 0
Q_vert = 0
q_down = 0
Q_down = 0

q_vert, Q_vert = vertical_wall_convection(Ts, Ti, t_open_t, l_vert, area_vert)
q_down, Q_down = free_convection_down(Ts, Ti, t_open_t, l_flat, area_flat)
q_total_open = q_vert + q_down
Q_total_open = Q_vert + Q_down

%% Forced Convection
t_forced_t= 45; % naturally convect same 
q_forced = 0
Q_forced = 0
q_forced, Q_forced = forced_convection(Ts, Ti, t_forced_t, .05 ,.025)


function[q_vert, Q_total ] = vertical_wall_convection(Ts,Ti,time, l_vert, area_vert)
global alp kin_v Pr k g B
% Vertical walls

%L_vert = .5 %meter
%l_vert = [.5 .25]
%area_vert = [ .25, .125] %meter2
Ra_vert = (g * B * (Ts - Ti).*l_vert.^3)/(alp * kin_v)
nusselt_vert = (.825 + (.387.*Ra_vert.^(1/6))/(1 + (.492/Pr)^(9/16))^(8/27)).^2

h_vert = nusselt_vert * k./l_vert
q_vert = h_vert .* area_vert * (Ts - Ti) % watts 
% copy code, for diff length
Q_total = q_vert * time 

end

function [q_forced, Q_total] = forced_convection(Ts, Ti, time, L_forced, area_forced)
velo2=200 %m/s
global kin_v Pr k 
%L_forced = .05 %char len
%area_forced = .025

% Recalculate kin_V, Pr
Rex = velo2 * L_forced/kin_v

nusselt_forced = .0296 * Rex ^.8 * Pr ^(1/3)
h_forced = nusselt_forced * k./L_forced
q_forced = h_forced .* area_forced * (Ts - Ti) % watts 
Q_total = q_forced * time
end


function [q_down, Q_total] = free_convection_down(Ts, Ti, time, L, area_flat)
%Facing down
%same L and A as bottom  plate
global alp kin_v k g B

Ra_down = (g * B * (Ts - Ti).*L.^3)/(alp * kin_v)
nusselt_down = .52 * Ra_down.^(1/5)
h_down = nusselt_down * k./L
q_down = h_down .* area_flat * (Ts - Ti) % watts 
Q_total = q_down * time
end


function [q_rad, Q_total] = radation(Ts, Ti, time, area_vert)
emiss= .5;
boltzman= 5.67e-8 %W/(m^2K^4)
q_rad= emiss*boltzman*(Ts^4-T1^4)*area_vert
end

function [q_free_up, Q_total] = free_convection_up(Ts, Ti, time, l_flat, area_flat)
% Flt plate facing up
global g B alp kin_v k
l_flat = [.75] % todo: L should be As/perimeter
area_flat= [.625]
Ra_flat = (g * B * (Ts - Ti).*l_flat.^3)/(alp * kin_v)
nusselt_flat = .15 * Ra_flat.^(1/3)
h_flat = nusselt_flat * k./l_flat
q_free_up = h_flat .* area_flat * (Ts - Ti) % watts 
Q_total = q_free_up * time

end

