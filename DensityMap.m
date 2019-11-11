clear all
Pressure = 1.11e+3 % kilapascals
Pressure = 6e3
Temp = -20 + 273
% 998 is gas
% if Q < .5 liquid
% pressure from 4 mpa to 60 mpa
% from 0 c to 250c 


p = 8e3 : .5e3 : 20e3; 

%p = 4e3 : 1e3 : 60e3 
t = 0 + 273:1:250 + 273;
h = 300e3 : 10e3: 700e3;

%[X,Y] = meshgrid(p,t)


for i = 1:length(p)
    for j = 1:length(h)
         D(i,j) = refpropm('D','H',h(j), 'P',p(i), 'CO2');
%          X(i,j) = refpropm('X','T',t(j), 'P',p(i), 'CO2');
         Q(i,j) = refpropm('Q','H',h(j), 'P',p(i), 'CO2');
%         X(i,j) = refpropm('X','H',h(j), 'P',p(i), 'CO2');
       
    end
end
%%
figure

% xlabel('Temperature (Kelvin)')
xlabel('Enthalpy')
ylabel('Pressure (kPa)')
title('CO2 Vapor Dome')
legend('Liquid', 'Gas', 'Super Critical', 'Mixed')

surf(h, p , D)
view(2)
xlabel('Enthalpy')
ylabel('Pressure (kPa)')

colorbar
caxis([0 400])
title('Density Color Map, Super Critical Capped at 400')
