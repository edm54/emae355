function [x] = plotVaporDome()
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Pressure = 1.11e+3 % kilapascals
Pressure = 6e3
Temp = -20 + 273
% 998 is gas
% if Q < .5 liquid
% pressure from 4 mpa to 60 mpa
% from 0 c to 250c 


p = 1e3 : .25e3 : 27.75e3; 

%p = 4e3 : 1e3 : 60e3 
t = 0 + 273:1:250 + 273;
h = 100e3 : 2e3: 700e3;

%[X,Y] = meshgrid(p,t)

gas1 = []
gas2 = []
liquid1 = [];
liquid2 = [];
sup1 = [];
sup2 = [];
mixed1 = [];
mixed2 = [];

for i = 1:length(p)
    for j = 1:length(h)
%         Q(i,j) = refpropm('Q','T',t(j), 'P',p(i), 'CO2');
%         X(i,j) = refpropm('X','T',t(j), 'P',p(i), 'CO2');
        Q(i,j) = refpropm('Q','H',h(j), 'P',p(i), 'CO2');
        X(i,j) = refpropm('X','H',h(j), 'P',p(i), 'CO2');
        if Q(i,j) < .01 
            liquid1 = [p(i)/1e3 liquid1 ];
            liquid2 = [h(j) liquid2 ];
        elseif Q(i,j) >= .01 && Q(i,j) <=1
            mixed1 = [p(i)/1e3 mixed1];
            mixed2 = [h(j) mixed2 ];
        elseif Q(i,j) == 998 || Q(i,j) == -998 || (Q(i,j) > 1 && Q(i,j)< 5)
            gas1 = [p(i)/1e3 gas1 ];
            gas2 = [h(j) gas2 ];
        else
            sup1 = [p(i)/1e3 sup1 ];
            sup2 = [h(j) sup2 ];
        end
    end
end
%%
x = figure
scatter(liquid2, liquid1  ,[], 'red', 'DisplayName', 'Liquid')
hold on
scatter(gas2, gas1, [], 'blue', 'DisplayName', 'Gas')
hold on
scatter(sup2, sup1, [], [.9 .5 0], 'DisplayName', 'Super Critical')
hold on
scatter(mixed2, mixed1, [], [.5 .5 .5], 'DisplayName', 'Super Mixed')
% xlabel('Temperature (Kelvin)')
xlabel('Enthalpy')
ylabel('Pressure (MPa)')
title('CO2 Vapor Dome')
legend('Liquid', 'Gas', 'Super Critical', 'Mixed')






end

