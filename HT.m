function [Tfluid2, enthalpy2] = HT(mdot,Pfluid1,Tfluid1, step_size, direction, Pfluid2, height)
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
Tco2 = Tfluid1; %initial CO2 temperature THIS INITIAL GUESS COMES FROM COMPRESSOR
Tgrad = 0.03562 ; % temperature gradient in K/m

%%
%will be used in convection calculation
% step_size = .5
% i = 1;
% Pfluid = linspace(60*1000,30*1000,6401);

if direction == 1 %This is for starting at top and going to bottom
%     for x = 0:step_size:L
%         if x == 0
%             Pfluid(i) = pressure;
%             Tsurf(i) = Ttop;
%             Tfluid(i) = Tco2;
%             i = i+1;
%         else
            Tground = Ttop + Tgrad*height;
            
            kfluid = refpropm('L','T',Tfluid1,'P',Pfluid1/1e3,'CO2'); %thermal conductivity of fluid
            rho = refpropm('D','T',Tfluid1,'P',Pfluid1/1e3,'CO2'); %density
            mu = refpropm('V','T',Tfluid1,'P',Pfluid1/1e3,'CO2'); %fluid viscosity Pa*s
            mus = refpropm('V','T',Tground,'P',Pfluid1/1e3,'CO2'); %surface temp visc
            
            Pr = refpropm('^','T',Tfluid1,'P',Pfluid1/1e3,'CO2');
            V = mdot/(rho*A);
            Re = rho*V*step_size/mu;
            Nu = (0.27*Re^(4/5))*Pr^1/3*(mu/mus)^.14;           
            h = Nu*kfluid/(2*r0_in);
            
            Rtot = (log(r2/r0_in)/k_steel)+(log(r3/r2)/k_conc)+...
                (log(r0_out/r3)/k_steel)+(1/(h*A));            
            
            
            enthalpy1 = refpropm('H','T',Tfluid1,'P',Pfluid1/1e3,'CO2');
            enthalpy2 = (Tground-Tfluid1)/(Rtot*mdot) + enthalpy1;
            
            %Pfluid2 = pressure_drop_down_one(mdot, Pfluid1, step_size, Tfluid1);
            
            Tfluid2 = refpropm('T','H',enthalpy2,'P',Pfluid2/1e3,'CO2');
            
%         end
%     end
elseif direction == 2 %Going from bottom up
%         for x = 0:.5:L
%             if x == 0
%                 Pfluid(i) = pressure;
%                 Tsurf(i) = Tbot;
%                 Tfluid(i) = Tco2;
%                 i = i+1;
%             else
             Tground = Tbot - Tgrad*(L-height);
                
                
            kfluid = refpropm('L','T',Tfluid1,'P',Pfluid1/1e3,'CO2'); %thermal conductivity of fluid
            rho = refpropm('D','T',Tfluid1,'P',Pfluid1/1e3,'CO2'); %density
            mu = refpropm('V','T',Tfluid1,'P',Pfluid1/1e3,'CO2'); %fluid viscosity Pa*s
            mus = refpropm('V','T',Tground,'P',Pfluid1/1e3,'CO2'); %surface temp visc
            
            Pr = refpropm('^','T',Tfluid1,'P',Pfluid1/1e3,'CO2');
            V = mdot/(rho*A);
            Re = rho*V*step_size/mu;
            Nu = (0.27*Re^(4/5))*Pr^1/3*(mu/mus)^.14;           
            h = Nu*kfluid/(2*r0_in);
            
            Rtot = (log(r2/r0_in)/k_steel)+(log(r3/r2)/k_conc)+...
                (log(r0_out/r3)/k_steel)+(1/(h*A)); 

                enthalpy1 = refpropm('H','T',Tfluid1,'P',Pfluid2/1e3,'CO2');
                enthalpy2 = (Tground-Tfluid1)/(Rtot*mdot) + enthalpy1;
                
                %[Pfluid2, pressure_loss2, gravity_gain2] = pressure_drop_up(mdot, Pfluid1);
            
                Tfluid2 = refpropm('T','H',enthalpy2,'P',Pfluid2/1e3,'CO2');
                
                
%             end
%         end
end

end
