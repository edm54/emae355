% EMAE 360 Real Otto/Atkinson Cycle Analysis Tool

%The purpose of this program is to set up a framework to solve for state
%conditions in a supercritical CO2 power cycle


%Author: Carter Waligura

%Date: 11/11/2019

clc
clear
close all

set(0, 'DefaultAxesFontWeight', 'normal', ...
      'DefaultAxesFontSize', 18, ...
      'DefaultAxesFontAngle', 'normal', ... 
      'DefaultAxesFontWeight', 'normal', ... 
      'DefaultAxesTitleFontWeight', 'bold', ...
      'DefaultAxesTitleFontSizeMultiplier', 1.2) ;
set(groot,'defaultLineLineWidth',3)


%% States

%1= input 
%2= after compressor
%3= after recuperator
%4= left bottom of heat exchanger hole
%5= after flowing through sand at bottom of hole
%6= right before turbine
%7= after turbine
%8= after recuperator part 2
%9= after cooler (same conditions as one with different mass flow)

% Trading Variables
P6= 20; % MPa
m2= 10; %kg/s

Eff_comp= .9;
Eff_turb= .9;
UA_re= 30; %might have to solve for these?
UA_cool= 30; 

% Knowns
P1= 6; %MPa
T1= 15+273.15; %K
m1= 3.24; %kg/s

m3= m1+m2;
P7=P1;


%Equations
%h1 = refpropm('H','T',T1, 'P',P1, 'CO2');
h1= h_T(T1,P2);

%s1= refpropm('S','T',T1, 'P',P1, 'CO2');
s1= s(T1,P1);

Pratio= P6/P1;

P5= f(friction, gravity, P6);
P4= f(darcy);
P3= f(friction, gravity, P4);
%Same with T's
T6=f(friction, T3);

P2= P3;

%h2_i= refpropm('H','S',s1, 'P',P1, 'CO2');
h2_i= hs(s1,P1);

W12= m3*(h2_i-h1);
h2= h1+ W12*Eff_comp/m3;

%T2= refpropm('T','H',h2, 'P',P2, 'CO2');
T2= T(h2,P2);

%h6= refpropm('H','T',T6, 'P',P5, 'CO2');
h6= h_T(T6,P5);
%Unknowns
% m3*(h2-h3)-m2*(h8-h7)=0
% m3*(h2-h3)- UA_re*((T2-T8)-(T3-T7))/ln((T2-T8)/(T3-T7))=0;
% m2*(h8-h1)- UA_cool*(T8-T1)=0;

%  T3=?;
%  T8=?;
% From unknowns
h3= h_T(T3,P2);
%h3= refpropm('H','T',T3, 'P',P2, 'CO2');
%s6= refpropm('S','T',T6, 'P',P6, 'CO2');
s6= s(T6, P6);

%h7_i= refpropm('H','S',s6, 'P',P7, 'CO2');
h7_i= h_s(s6,P7);

W67= m2*(h7_i-h6);
h7= h6+W67*Eff_turb/m2;
T7= refpropm('T','H',h7, 'P',P7, 'CO2');


function h_s = h_s(s, P)
    h_s = refpropm('H','S',s, 'P',P, 'CO2');
    %return
end

function h_T = h_T(T, P)
    h_T = refpropm('H','T',T, 'P',P, 'CO2');
    %return
end

function s = s(T, P)
    s = refpropm('H','T',T, 'P',P, 'CO2');
    %return
end

function T=T(h,P)
    T= refpropm('T','H',h, 'P',P, 'CO2');
    %return
end

