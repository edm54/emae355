clear all
clc
close all

% Set global variables for access within functions
% Generally not best practice but will work here
global alp kin_v  Pr k g B Dynamic_V rho
g = 9.81; %gravity, m/s 

Ts = 160 + 273;
Ti = 25 + 273 ;%77 f, room temperature
T_average = (Ts + Ti)/2;
B = 1/T_average;

% Prandtl number
Pr = refpropm('^','T',Ti,'P',101, 'Air.ppf');

% Dynamic viscosity [Pa*s], Mu
Dynamic_V = refpropm('V','T',Ti,'P',101, 'Air.ppf');

% Kinematic viscosity [cm^2/s], Nu
%Kin_v = Dynamic_v/rho
kin_v = refpropm('$','T',Ti,'P',101, 'Air.ppf');
kin_v = kin_v /100^2; % Convert to square meters

% Thermal diffusivity [cm^2/s], alpha
alp = refpropm('%','T',Ti,'P',101, 'Air.ppf');
alp = alp/ 100^2; % Convert to square meters

% Thermal conductivity [W/(m K)]
k = refpropm('L','T',Ti,'P',101, 'Air.ppf');

% Density
rho = refpropm('D','T',Ti,'P',101, 'Air.ppf');

faces = 4 ;
t_total = 7 * 60; % second

%% Inout Paramaters

area_flat= .891; % m^2 for ONE FLAT Face
perimeter_flat = 3.82; %meter
% Characteristic Length
l_flat = area_flat / perimeter_flat;

% Total area of 4 vertical walls
area_vert = .65/2; %meter^2;
perimeter_vert = 8.32/2;
% Characteristic length from drawings given
l_vert = .085; %m

%% Free Convection From Vertical Walls
t_vert = 7*60; %sec
[q_vertical, Q_vertical] = vertical_wall_convection(Ts, Ti, t_vert, l_vert,  area_vert);

%% Convection from Upper Surface of Hot Plate
t_bottom_convection= 3*60;
[q_bottom, Q_bottom] = free_convection_up(Ts, Ti, t_bottom_convection, l_flat, area_flat);

%% Convection while top plate is floating (i.e top and bottom surfaces fully exposed)
t_float_t= 30; % seconds

% Convection losses from top and bottom of the hot plate
[q_float_b, Q_float_b] =  free_convection_down(Ts, Ti, t_float_t, l_flat, area_flat);
[q_float_t, Q_float_t] = free_convection_up(Ts, Ti, t_float_t, l_flat, area_flat);

q_float = q_float_b  + q_float_t;
Q_float = Q_float_b + Q_float_t;

%% Convection from the bottom surface of a hot plate
t_open_t= 2.5*60; % seconds
[q_top, Q_top] = free_convection_down(Ts, Ti, t_open_t, l_flat, area_flat);

%% Forced convection due to spraying surface with compressed air
t_forced_t= 45; % naturally convect same 
[q_forced, Q_forced] = forced_convection(Ts, Ti, t_forced_t, .1*l_flat, .1*area_flat);

%% Radiation losses
t_rad = 7 * 60;
[q_rad, Q_rad] = radiation(Ts, Ti, t_rad, 2.*area_vert, .45);

%% Forced Convection from air flowing through platen
t_f = 5; %sec
[q_forced_pipe, Q_forced_pipe] = forced_pipe_convection(Ts, Ti,t_f);

%% New Calcs
t_rad_n = 6 * 60;
t_vert_n = 6 * 60;
t_bottom = 2 * 60;
t_top = 1.5 * 60;

% Solec LO/MIT 1 has emissivity .145
[q_rad_n, Q_rad_n] = radiation(Ts, Ti, t_rad_n , 2.*area_vert, .145);
[q_vertical_n, Q_vertical_n] = vertical_wall_convection(Ts, Ti, t_vert_n,  l_vert, area_vert);
[q_top_n, Q_top_n] = free_convection_down(Ts, Ti, t_top, l_flat, area_flat);
[q_bottom_n, Q_bottom_n] = free_convection_up(Ts, Ti, t_bottom, l_flat, area_flat);

%% Savings calculations
Q_old_total = Q_forced + Q_bottom + Q_float ...
                 + Q_top + 2*Q_vertical + Q_rad + Q_forced_pipe;   
             
Q_old_top = Q_float_t + Q_top + .5 * Q_rad + Q_forced + Q_vertical + Q_float_b;
Q_old_bot = Q_bottom + Q_vertical + .5 * Q_rad + Q_forced_pipe;
     
Q_new_total = Q_forced  + Q_bottom_n + Q_float ...
                 + Q_top_n + 2*Q_vertical_n + Q_rad_n;       
Q_new_top = Q_float_t + Q_top_n + .5 * Q_rad_n + Q_forced + Q_vertical_n +Q_float_b;
Q_new_bot = Q_bottom_n + Q_vertical_n + .5 * Q_rad_n;          

savings = (Q_old_total - Q_new_total)/Q_old_total     
Q_loss = Q_old_total - Q_new_total;
cycles = 24 * 60 /7;

heat_loss_daily_original = Q_old_total * cycles ;
heat_loss_daily_new = Q_new_total * cycles ;
heat_loss_daily_savings = heat_loss_daily_original - heat_loss_daily_new;


%%
total_toys_per_machine = 30 * 24 * 60 /7;
time_with_new_cycle = total_toys_per_machine/30 * 6/60;

% Hours of productivity saved in a day with solution
time_saved = 24 - time_with_new_cycle;

% Min number of machines to turn off one
min_machines = 24/time_saved;

daily_amount_saved_pm = (heat_loss_daily_original - heat_loss_daily_new)/7;
cycles = total_toys_per_machine/30;

%% Payback:
% an electricity rate of 6¢/kWhr for use in the payback period estimate
machines = 1; % Does not effect final answer
cost_pin = 4.72; 
surface_area_needed = 1.3 * machines; % would be .65 but assume some overspray/leakage
bucket_cost = 760; % Dollars per bucket 
coverage = 116.1288; %m2 this is 1250 sqft
coating_cost = bucket_cost * surface_area_needed/coverage;
possible_machines = coverage/1.3;

soln_cost = 30 * machines * cost_pin + coating_cost + 70; % Dollars (240 pins plus a bucket) for 8 machines
e_cost = .06; % Dollars/kWhr

% Convert Joules to kWhr
energy_save_day = heat_loss_daily_savings/3.600e6;
daily_money_saved = machines * .06 * energy_save_day;
days_to_even = soln_cost/daily_money_saved;
%%
full_day_old = 7 * cycles * Q_old_total;
new_cycles = 7 * cycles/ 6;
full_day_new = 6 * new_cycles * Q_new_total;
savings_new = (full_day_old - full_day_new)/full_day_old;


%% Plotting
figure
%labels = {'Forced', 'Bottom' , 'Float', 'Top', 'Vertical', 'Radiation', 'Forced Pipe'}
h1 = pie([Q_forced Q_bottom  Q_top  Q_float 2*Q_vertical Q_rad Q_forced_pipe]);
colormap([230/256 230/256 .98; 52/256 190/256 235/256; ...
          70/256 130/256 235/256; 0/256 50/256 190/256; 0 18/256 105/256;  0/256 0/256 55/256; 1 1 1]);
leg = legend('Forced: Spraying Air', ...
       'Bottom Plate: Top Surface' , 'Top Plate: Bottom Surface',  ...
       'Top Plate: Both Surfaces', 'Free Convection: Vertical Walls', ...
       'Radiation', 'Forced: Air Through Platten');
leg.FontSize = 14;   
set(h1(2:2:end),'FontSize',12);
title('Original Heat Transfer Distribution');
ax = gca;
ax.TitleFontSizeMultiplier = 1.6;

figure
pie([q_forced q_bottom q_float ...
                 q_top  2*q_vertical  q_rad, q_forced_pipe]);
legend('Forced: Spraying Air', 'Bottom Plate: Top Surface' , 'Top Plate: Both Surfaces', ...
        'Top Plate: Bottom Surface', 'Free Convection: Vertical Walls', 'Radiation', ...
        'Forced: Air Through Platten');
title('Original q Distribution');
colormap([230/256 230/256 .98; 52/256 190/256 235/256; ...
          70/256 130/256 235/256; 0/256 50/256 190/256; 0 18/256 115/256;  0/256 0/256 55/256; .97 .97 1]);

figure
h2 = pie([Q_forced Q_bottom_n  Q_top_n Q_float 2*Q_vertical_n  Q_rad_n Q_loss],...
    [0 0 0 0 0 0 1]);
colormap([230/256 230/256 .98; 52/256 190/256 235/256; ...
          70/256 130/256 235/256; 0/256 50/256 190/256; 0 18/256 110/256;  0/256 0/256 55/256; .97 .97 1])
leg = legend('Forced: Spraying Air', 'Bottom Plate: Top Surface' , ...
       'Top Plate: Bottom Surface', 'Top Plate: Both Surfaces', ...
       'Free Convection: Vertical Walls', 'Radiation', 'Savings');
leg.FontSize = 14;   
set(h2(2:2:end),'FontSize',14);
title('Heat Transfer Distribution With Solution');
ax = gca;
ax.TitleFontSizeMultiplier = 1.6;
%%
figure
h2 = pie([Q_forced Q_bottom_n  Q_top_n Q_float 2*Q_vertical_n  Q_rad_n ]);
colormap([230/256 230/256 .98; 52/256 190/256 235/256; ...
          70/256 130/256 235/256; 0/256 50/256 190/256; 0 18/256 110/256;  0/256 0/256 55/256])
leg = legend('Forced: Spraying Air', 'Bottom Plate: Top Surface' , ...
       'Top Plate: Bottom Surface', 'Top Plate: Both Surfaces', ...
       'Free Convection: Vertical Walls', 'Radiation');
title('Heat Transfer with Solution, No Savings')

%%
% Free convection from vertical walls
function[q_vert, Q_total ] = vertical_wall_convection(Ts,Ti,time, l_vert, area_vert)
global alp kin_v Pr k g B
    Ra_vert = (g * B * (Ts - Ti).*l_vert.^3)/(alp * kin_v);
    nusselt_vert = (.825 + (.387.*Ra_vert.^(1/6))/(1 + (.492/Pr)^(9/16))^(8/27)).^2;
    h_vert = nusselt_vert * k./l_vert;
    q_vert = h_vert .* area_vert * (Ts - Ti); % watts 
    Q_total = q_vert * time ;
end

% Forced convection from spraying air over surface
function [q_forced, Q_total] = forced_convection(Ts, Ti, time, L_forced, area_forced)
    global kin_v Pr k 
    velo2=25; %m/s
    Rex = velo2 * L_forced/kin_v;
    nusselt_forced = .0296 * Rex ^.8 * Pr ^(1/3);
    h_forced = nusselt_forced * k./L_forced;
    q_forced = h_forced .* area_forced * (Ts - Ti);
    Q_total = q_forced * time;
end

% Convection from the bottom surface of a hot plate
function [q_down, Q_total] = free_convection_down(Ts, Ti, time, L, area_flat)
    global alp kin_v k g B
    Ra_down = (g * B * (Ts - Ti).*L.^3)/(alp * kin_v);
    nusselt_down = .52 * Ra_down^(1/5);
    h_down = nusselt_down * k./L;
    q_down = h_down .* area_flat * (Ts - Ti); % watts 
    Q_total = q_down * time;
end


function [q_rad, Q_total] = radiation(Ts, Ti, time, area_vert, emiss)
    boltzman= 5.67e-8; %W/(m^2K^4);
    q_rad= emiss*boltzman*(Ts^4-Ti^4)*area_vert;
    Q_total  = q_rad * time;
end

% Convection from the top surface of a hot plate
function [q_free_up, Q_total] = free_convection_up(Ts, Ti, time, l_flat, area_flat)
    global g B alp kin_v k
    Ra_flat = (g * B * (Ts - Ti).*l_flat.^3)/(alp * kin_v);
    nusselt_flat = .15 * Ra_flat.^(1/3);
    h_up = nusselt_flat * k./l_flat;
    q_free_up = h_up .* area_flat * (Ts - Ti); % watts 
    Q_total = q_free_up * time;
end

function [q_pipe_convection, Q_pipe_convection] = forced_pipe_convection(Ts, Ti, time)
    global rho Pr k
    Dynamic_V = refpropm('V','T',Ti,'P',101, 'Air.ppf');
    V = 25; %m/s, average
    D = 5/1000; %m --> 5mm
    length = 30/1000; %m --> 30 mm
    Re = rho * V * D/Dynamic_V;

    Nu = .023 * Re^.8 * Pr^.4;
    h = Nu * k/ D; % w/m^2 K
    area_affected = pi * D * length;
    q_pipe_convection = 30 *h .* area_affected * (Ts - Ti); % watts 
    Q_pipe_convection = q_pipe_convection * time;
end