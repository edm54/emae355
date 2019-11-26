% EMAE 360 Real Otto/Atkinson Cycle Analysis Tool

%The purpose of this program is to set up a framework to solve for state
%conditions in a supercritical CO2 power cycle

%This script does not contain a recuperator in the analysis

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
%3= left bottom of heat exchanger hole
%4= after flowing through sand at bottom of hole
%5= right before turbine
%6= after turbine
%7= after cooler (same conditions as one with different mass flow)

%Kept same variable naming convention as recuperator for clarity


%Energies in J
%Pressures inn kPa
%Temperatures in K



%Initial Guesses: Will need to look at chart
Eff_comp= .8;
Eff_turb= .8;
%UA_cool= 30; 

% Knowns
P1= 6*1e6; %kPa
T1= 15+273.15; %K
m1= 3.24; %kg/s

h1= h_T(T1,P1); 

P7=P1;

s1= s(T1,P1);

% Trading Variables
P2= (23:3:29)*1e6; % Pa
m2= 7:.2:14; %kg/s
Pratiocomp= P2/P1;

figure
legendstr= {};

T2=ones(length(P2), length(m2))*303; %initial T2 estimate
for i= 1:length(P2)
    for j=1:length(m2)
          m3(j)= m1+m2(j);
        %Big Function
%               [P5(i,j), T5(i,j), P3(i,j), T3(i,j), P4(i,j), T4(i,j)]...
%             = findP5(P2(i),m3(j));
        try
            if j==1
            [P5(i,j), T5(i,j), P3(i,j), T3(i,j), P4(i,j), T4(i,j)]...
                = findP5(P2(i),m3(j), T2(i,j));
            else
                
            [P5(i,j), T5(i,j), P3(i,j), T3(i,j), P4(i,j), T4(i,j)]...
                = findP5(P2(i),m3(j), T2(i,j-1));
            end
        catch ME
%           P5(i,j)=NaN;
%           T5(i,j)=NaN;
          P5(i,j)=P5(i,j-1);
          T5(i,j)=T5(i,j-1);
        end
        
%         if j>1
%             if Pturb(i,j-1)<0
%                 k=j-1;
%                 [P5(i,j), T5(i,j)] = findP5(P2(i),m3(k));
%             end
%         end
                
        h5(i,j)= h_T(T5(i,j),P5(i,j));
        s5(i,j)= s(T5(i,j), P5(i,j));

        % P3= f(friction, gravity,heat, P2);
        % P4= f(darcy);
        % P5= f(friction, gravity, P4);
        % %Same with T's
        % T5=f(friction, T3);

        %Need P6 and T6


        h2_i(i)= h_s(s1,P2(i));

        W12(i,j)= m3(j)*(h2_i(i)-h1);
        h2(i,j)= h1+ W12(i,j)*Eff_comp/m3(j);

        %T2= refpropm('T','H',h2, 'P',P2, 'CO2');
        T2(i,j)= T(h2(i,j),P2(i));

        %h6= refpropm('H','T',T6, 'P',P5, 'CO2');
        %Unknowns
        % m2*(h8-h1)- UA_cool*(T8-T1)=0;

        P6=P1; %=P7

        % From unknowns
        %h3= h_T(T3,P2);
        %h3= refpropm('H','T',T3, 'P',P2, 'CO2');
        %s6= refpropm('S','T',T6, 'P',P6, 'CO2');

        %h7_i= refpropm('H','S',s6, 'P',P7, 'CO2');
        h6_i(i,j)= h_s(s5(i,j),P6);


        W56(i,j)= m2(j)*(h6_i(i,j)-h5(i,j));
        h6(i,j)= h5(i,j)+W56(i,j)*Eff_turb/m2(j);
        T6(i,j)= T(h6(i,j), P6);
        UA_cool(i,j)= m2(j)*(h6(i,j)-h1)/(T6(i,j)-T1);

        Pturb(i,j)= m2(j)*(h5(i,j)-h6(i,j));
        Pcomp(i,j)= m3(j)*(h2(i,j)-h1); %should be in watts
       % Pdiff(i,j)=Pturb(i,j)-Pcomp(i,j);
    end
    disp(P2(i)/1e6)
    Pdiff(i,:)= Pturb(i,:)-Pcomp(i,:);
    plot(m2,Pdiff(i,:)/1000)
    hold on
    legendstr{i}=strcat('P2=', num2str(P2(i)/1000));
    %legendstr= {num2str(P2(i)/1000) legendstr};
end
legend(legendstr)
xlabel('Mass Flow Through Turbine')
ylabel('Net Power (kW)')
title('Power vs Mass Flow with varying Initial Pressures')
%axis([ m2(1) m2(end) 0 100])
grid

figure
for w=1:length(P2)
plot(m2, Pturb(w,:), '--')
%legend(legendstr)
hold on
plot(m2, Pcomp(w,:), '-')
end
legend('turb=20',  'comp=20', 'turb=40','comp=40', 'turb=60', 'comp=60')
grid
%legend(legendstr)

%% Choosing Machinery
max= max(max(Pdiff));
[I,J]=find(Pdiff==max);

%Find machinery stuff using the best pressure and mdot case

Tavg_comp= (T1+T2(I,J))/2;
Pavg_comp= (P1+P2(I))/2;

Tavg_turb= (T5(I,J)+T6(I,J))/2;
Pavg_turb= (P5(I,J)+P6)/2;

%Isentropic Head
% His_comp= gammas(Tavg_comp,Pavg_comp)/(gammas(Tavg_comp,Pavg_comp)-1)*Pavg_comp...
%     /rho(Tavg_comp, Pavg_comp)*((P2(I)/P1)^((gammas(Tavg_comp,Pavg_comp)-1)/...
%     gammas(Tavg_comp,Pavg_comp))-1); %equation from textbook
% 
% His_turb= -gammas(Tavg_turb,Pavg_turb)/(gammas(Tavg_turb,Pavg_turb)-1)*Pavg_turb...
%     /rho(Tavg_turb, Pavg_turb)*((P7/P5(I,J))^((gammas(Tavg_turb,Pavg_turb)-1)/...
%     gammas(Tavg_turb,Pavg_turb))-1);
His_comp= (h2(I,J)-h1)/9.81;
His_turb= (h5(I,J)-h6(I,J))/9.81;



Vdot_comp= m3(J)*rho(T1,P1);
Vdot_turb= m2(J)*rho(T6(I,J),P7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 1
N1=3600;

Ns1comp= N1*sqrt(Vdot_comp)/His_comp^.75;
Ns1turb= N1*sqrt(Vdot_turb)/His_turb^.75;
Ns1turb2=N1*2*pi/60*sqrt(Vdot_turb)/(His_turb*9.81)^.75; %WHICH ONE IS RIGHT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Case 2
%Nstarget= 300; 
N2=1000:100:30000; %set range of RPMs to try and find Ns?
D2= .05:.05:1;

Ns2comp= N2*sqrt(Vdot_comp)/His_comp^.75;
Ns2turb= N2*sqrt(Vdot_turb)/His_turb^.75;

Ds2comp= D2*His_comp^.25/sqrt(Vdot_comp);
Ds2turb= D2*His_turb^.25/sqrt(Vdot_turb);


figure
plot(N2, Ns2comp)
hold on
plot(N2, Ns2turb)
xlabel('RPM (N)')
ylabel('Specific RPM (Ns)')
legend('specN Pump', 'specN Turb')
title(['Specific RPM vs RPM at P= ' num2str(P2(I)/1e6) 'MPa and m=' num2str(m2(J)) 'kg/s'])

figure
plot(D2, Ds2comp)
hold on
plot(D2, Ds2turb)
xlabel('Diameter (m)')
ylabel('Specific Diameter (Ds)')
legend('specD Pump', 'specD Turb')
title(['Specific Diameter vs Diameter at P= ' num2str(P2(I)/1e6) 'MPa and m=' num2str(m2(J)) 'kg/s'])


f1=figure;
pl1=plot([100, 310, 40 100] , [1 .85 2 1], 'r-'); %Range for .8 eff
hold on
pl2=plot([20, 100, 800 20] , [4 .8 .8 4], 'b-'); %Range for .7 eff
pl123=plot([100 2000 6 100]  ,[.6 .6 7.5 .6], 'g-'); %range for .6 eff
% plot(meshgrid(Ns2turb, Ds2turb))
for a=1:length(Ns2turb)
    for b=1:length(Ds2turb)
    plot(Ns2turb(a), Ds2turb(b), 'xk')
    hold on
    end
end
pl1=plot([100, 310, 40 100] , [1 .85 2 1], 'r-'); %Range for .8 eff
hold on
pl2=plot([20, 100, 800 20] , [4 .8 .8 4], 'b-'); %Range for .7 eff
pl123=plot([100 2000 6 100]  ,[.6 .6 7.5 .6], 'g-'); %range for .6 eff
hold on
% copyobj(pl1, f1)
% copyobj(pl2, f1)
set(gca,'XScale','log','YScale','log')
title('Turbine Trade Study')
xlabel('Specific Speed (Ns)')
ylabel('Specific Diameter (Ds)')
legend('.8 eff', '.7 eff', '.6 eff',  'Calc')

f2=figure;
pl3=plot([60, 175, 350, 60] , [2.85 .85 .8 2.85], '-r'); %range for .8 for axial 
hold on
pl4=plot([45, 150, 800, 45] , [4 .8 .6 4], '-b'); %range for .7 
pl15= plot([2000 30 150 2000]   , [.5 5 .65 .5], 'g'); %range for .6
for a=1:length(Ns2comp)
    for b=1:length(Ds2comp)
    plot(Ns2comp(a), Ds2comp(b), 'xk')
    hold on
    end
end
pl3=plot([60, 175, 350, 60] , [2.85 .85 .8 2.85], '-r'); %range for .8 for axial 
hold on
pl4=plot([45, 150, 800, 45] , [4 .8 .6 4], '-b'); %range for .7 
pl15= plot([2000 30 150 2000]   , [.5 5 .65 .5], 'g'); %range for .6% 

title('Compressor Trade Study')
xlabel('Specific Speed (Ns)')
ylabel('Specific Diameter (Ds)')
legend('.8 eff', '.7 eff', '.6 eff', 'Calc')
set(gca,'XScale','log','YScale','log')
%Maybe multiple by 129

%% P h diagram

hvec= [h1 h2(I,J) h_T(T3(I,J),P3(I,J)) h_T(T4(I,J), P4(I,J)) h5(I,J) h6(I,J) h1]/1e3;
Pvec= [P1 P2(I) P3(I,J) P4(I,J) P5(I,J) P6 P1];

%figure
fig= plotVaporDome();
hold on
for g=1:6
    plot([hvec(g) hvec(g+1)], [Pvec(g)/1e6 Pvec(g+1)/1e6])
    hold on
end
legend('Liquid' ,'Gas', 'Supercritical' ,'Mixed' ,'1-2', '2-3', '3-4', '4-5' ,'5-6', '6-1')
%plot(hvec/1000, Pvec/1e6)
xlabel('Enthaply (kJ/kg)')
ylabel('Pressure (MPa)')
title(['P-h Power Cycle Diagram at P2= ' num2str(P2(I)/1e6) 'MPa and m2=' num2str(m2(J)) 'kg/s'])

%% Functions

function h_s = h_s(s, P) %J/kg
    h_s = refpropm('H','P',P/1e3,'S',s, 'CO2');
    %return
end

function h_T = h_T(T, P) %J/kg
    h_T = refpropm('H','T',T, 'P',P/1e3, 'CO2');
    %return
end

function s = s(T, P) %J/kgK
    s = refpropm('S','T',T, 'P',P/1e3, 'CO2');
    %return
end

function T=T(h,P) %K
    T= refpropm('T','H',h, 'P',P/1e3, 'CO2');
    %return
end

function gammas= gammas(T,P) 
    gammas= refpropm('K','T', T, 'P', P/1e3, 'CO2');
end

function rho= rho(T,P) %kg/m^3
    rho= refpropm('D','T', T, 'P', P/1e3, 'CO2');
end
