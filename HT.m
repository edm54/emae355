function Tfluid = HT(mdot,pressure,temp,tob)
% The function takes a mass flow rate, initial fluid pressure, and
% temperature coming out of the pump. 
% tob is a factor for going up or down the pipe
% tob = 1 goes from the compressor/pump down the well and tob = 2 goes up
% the pipe
% The result is the fluid temperature at each specified depth, accounting
% for changes in pressure. 

%% Values
r0_in = 0.0610; %m
r0_out = 0.1397; %m
t_steel = 0.0159; %thickness assumed for at least 60MPa pressure rated piping (McMaster Carr)
t_conc = 0.0628;
k_conc = 1.4; % W/mK @300K
k_steel = 60.5; %W/mK @300K in textbook tables
r2 = (r0_in+t_steel);
r3 = r0_in+t_steel+t_conc;
L = 3200; %depth 
R = 188.92; %specific gas constant CO2
A = (r0_in)^2*pi/2; %pipe area

Tbot = 129+273; % initial pipe surface temperature (ground temp) at bottom
Ttop = 15+273; %top temp
Tco2 = temp; %initial CO2 temperature
Tgrad = 0.30875 ; % temperature gradient in K/m

%%
step_size = .5
i = 1;
Pfluid = linspace(60*1000,30*1000,6401);

if tob == 1
    for x = 0:step_size:L
        if x == 0
            Pfluid(i) = pressure;
            Tsurf(i) = Ttop;
            Tfluid(i) = Tco2;
            i = i+1;
        else
            Tsurf(i) = Ttop + Tgrad*x;
            
            Rtot = (log(r2/r0_in)/k_steel)+(log(r3/r2)/k_conc)+(log(r0_out/r3)/k_steel);
            
            enthalpy(i-1) = refpropm('H','T',Tfluid(i-1),'P',Pfluid(i),'CO2');
            enthalpy(i) = (Tsurf(i)-Tfluid(i-1))/(Rtot*mdot) + enthalpy(i-1);
            
            Pfluid(i) = pressure_drop_down_one(mdot, Pfluid(i-1), step_size, Tfluid(i-1));
            
            Tfluid(i) = refpropm('T','H',enthalpy(i),'P',Pfluid(i),'CO2');
            
            i = i+1;
        end
    end
else if tob == 2
        for x = 0:.5:L
            if x == 0
                Pfluid(i) = pressure;
                Tsurf(i) = Tbot;
                Tfluid(i) = Tco2;
                i = i+1;
            else
                Tsurf(i) = Tbot - Tgrad*x;
                
                Rtot = (log(r2/r0_in)/k_steel)+(log(r3/r2)/k_conc)+(log(r0_out/r3)/k_steel);

                enthalpy(i-1) = refpropm('H','T',Tfluid(i-1),'P',Pfluid(i),'CO2');
                enthalpy(i) = (Tsurf(i)-Tfluid(i-1))/(Rtot*mdot) + enthalpy(i-1);
                
                [Pfluid(i), pressure_loss(i), gravity_gain(i)] = pressure_drop_up(mdot, Pfluid(i-1));
            
                Tfluid(i) = refpropm('T','H',enthalpy(i),'P',Pfluid(i),'CO2');
                
                i = i+1;
            end
        end
    end
end
