function Tfluid = HT(mdot,pressure,temp_b)
% The function takes a mass flow rate, initial fluid pressure, and
% temperature.
% The result is the net heat transfer.

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

%% Boundary Conditions
T0 = 288;
Tco2 = temp_b; %initial CO2 temperature, in the "slice" below
P0 = 6;
Tf = 402;
Pf = 32.5;
Tm = (15 + 129)/2 + 273; %mean temperature 
Tgrad = 0.30875 ; % temperature gradient in K/m

%%
i = 1;
if i = 1:
    e(i) = 230; % random assumed value for enthalpy
    T(i) = 288; % assumed value for fluid temperature 

    enthalpy(i) = refpropm('H','T',Tco2,'P',pressure,'CO2');
    Pfluid(i) = pressure;
    Tsurf(i) = T0;
    Tfluid(i) = Tco2;
    i = i+1;
    
else
    Tsurf(i) = T0 + Tgrad*x;
    Pfluid(i) = pressure_drop_down(mdot, pressure);

    % ITERATION NEEDED FOR TEMPERATURE OF FLUID
    while diff > 0 
        e(i) = (T(i) - Tfluid(i-1))/mdot + enthalpy(i-1); % steady state constant heat transfer rate
        enthalpy(i) = refpropm('H','T',T(i),'P',Pfluid(i),'CO2');
        diff = e(i)- enthalpy(i);
    end

    Tfluid(i) = refpropm('T','H',enthalpy(i),'P',Pfluid(i),'CO2');
    q(i) = (2*pi*x*(Tsurf(i)-Tfluid(i)))/((log(r2/r0_in)/k_steel)+(log(r3/r2)/k_conc)+(log(r0_out/r3)/k_steel));
    i = i+1;
    
end


end
