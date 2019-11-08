clear all
Pressure = 1.11e+3 % kilapascals
Pressure = 6e3
Temp = -20 + 273
% 998 is gas
% if Q < .5 liquid
% pressure from 4 mpa to 60 mpa
% from 0 c to 250c 


p = 4e3 : 1e2 : 60e3 

%p = 4e3 : 1e3 : 60e3 
t = 0 + 273:.5:250 + 273
%[X,Y] = meshgrid(p,t)

gas1 = []
gas2 = []
liquid1 = [];
liquid2 = [];
sup1 = [];
sup2 = [];


for i = 1:length(p)
    for j = 1:length(t)
        Q(i,j) = refpropm('Q','T',t(j), 'P',p(i), 'CO2');
        X(i,j) = refpropm('X','T',t(j), 'P',p(i), 'CO2');
        if Q(i,j) < .5
            liquid1 = [p(i) liquid1 ];
            liquid2 = [t(j) liquid2 ];
        elseif Q(i,j) == 998 || Q(i,j) == -998 || (Q(i,j) > .5 && Q(i,j)< 5)
            gas1 = [p(i) gas1 ];
            gas2 = [t(j) gas2 ];
        else
            
            sup1 = [p(i) sup1 ];
            sup2 = [t(j) sup2 ];
        end
    end
end
%%
scatter(liquid2, liquid1  ,[], 'red', 'DisplayName', 'Liquid')
hold on
scatter(gas2, gas1, [], 'blue', 'DisplayName', 'Gas')
hold on
scatter(sup2, sup1, [], [.9 .5 0], 'DisplayName', 'Super Critical')
title('CO2 Vapor Dome')
xlabel('Temperature (Kelvin)')
ylabel('Pressure (kPa)')
legend


