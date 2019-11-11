%% Darcy
clear all
k = 1e-12 % 
A = 20 % m^2
L = 112 % m
Pa = 32.5e6 % MPa
Tb = 129  + 273 %k 
Dynamic_V = refpropm('V','T',Tb,'P', Pa/1e3, 'CO2');
density = refpropm('D','T',Tb,'P',Pa/ 1e3, 'CO2');
mass_flow = 3.24 % kg/s
%m_flow = 3.24 : .1 : 5.24

i = 1
for m_flow = 3.24 : .1 : 10
%for i = length(m_flow)

    Q(i) = m_flow / density
    delta_p = Q * Dynamic_V * L /(k * A)
    %Pb = Q * Dynamic_V * L /(k * A) + Pa
    %Pa -  Q * Dynamic_V * L /(k * A)
    i = i + 1
end

figure
plot(3.24 :.1: 10, delta_p)
xlabel('Mass Flow Rate (kg/s)')
ylabel('Pressure Difference')
title('Mass Flow Rate vs. Pressure Difference')
%%
i = 1
pb1 = 32.5e6 : .5e6 : 40e6
for Pb = 32.5e6 : .5e6 : 40e6
%for i = length(m_flow)

   
    Q1(i) = (k * A* ( Pb - Pa)/(Dynamic_V * L))* density
    %Pb = Q * Dynamic_V * L /(k * A) + Pa
    %Pa -  Q * Dynamic_V * L /(k * A)
    i = i + 1
end

figure
plot(pb1, Q1)
xlabel('Input Pressure (MPa)')
ylabel('Mass Flow Rate (kg/s)')
title('Mass Flow Rate vs. Input Pressure')

