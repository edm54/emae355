clc 
clear
close all

% Initial conditions
cond_name = ['1: T=-20 ', '2: T=-10 ', '3: T=0 ', '4: T=10 ', '5: T=20 ', '6: T=30 ' ...
    '7: T=40 ', '8: T=50 ', 'Startup: T=20 '];
Ambient= [ -20, -10, 0, 10, 20, 30, 40, 50, 20];
T1 = [ -20,-10, 0 ,10, 20, 30, 40, 50, 20];
T2 = [-14, -2, 10,23, 38, 57, 71, 85, 27];
T3 = [ 33,39,47, 59, 71,85,104, 122, 79];
T4 = [ 200,200,200,200,200,200,200,200, 200];
T5 = [81, 81,87, 96, 103, 113, 133, 150, 184];
T6 = [ -2, 3, 15, 28,43, 62, 77, 91, 33];
P1 = [ 3.3,3.3, 3.9, 5, 6.2, 7.4, 9.2, 10, 6.2];
P2 = [ 14, 15, 16, 18, 20, 20.7, 19.9, 17.8, 11] ;
m_dot = [ 4.9, 4.9, 5.2, 5.8, 6.4, 6.5, 5.5, 4, 3];

% Convert temperatures from Celcius to Kelvin
T1K = T1+ 273.15;
T2K = T2 + 273.15;
T3K = T3+ 273.15;
T4K = T4 + 273.15;
T5K = T5+ 273.15;
T6K = T6 + 273.15;
P1k = P1 * 1000;
P2k = P2 * 1000;

% Determine initial Enthalpy, Entropy
for state = 1:length(T1)
    H_a(state) = refpropm('H','T',T2K(state),'P',P2k(state),'CO2');
    S_a(state) = refpropm('S','T',T2K(state),'P',P2k(state),'CO2');
    
    H_b(state) = refpropm('H','T',T4K(state),'P',P2k(state),'CO2');
    S_b(state) = refpropm('S','T',T4K(state),'P',P2k(state),'CO2');  
end

% Iteration algorithm to find P_choke for each valve
for state = 1:length(T1)
    % Some of the states for T3 do not choke, so we use a try-catch 
    try
    % Guess intial choke pressure 
    Pchoke_a(state) = min(13000, P2k(state)-1000);% kPa
    Ma_a(state) = .5;
        while abs(Ma_a(state)-1) > .001
            SS_a(state) = refpropm('A','P',Pchoke_a(state),'S',S_a(state),'CO2');
            H2_a(state) = refpropm('H','P',Pchoke_a(state),'S',S_a(state),'CO2');            
            V2_a(state) = sqrt(2*(H_a(state)-H2_a(state)));
            Ma_a(state) = V2_a(state)/SS_a(state);
            
            % Increase P_Choke if Ma>1
            Pchoke_a(state) = Pchoke_a(state) + (Ma_a(state) - 1) * 50;
        end
    catch ME
        % Sets the choke pressure equal to the final pressure when the mach
        % state can not be reached. This is because the flow pressure can
        % now communicate up to the choke without the Mach 1 flow
        Pchoke_a(state) = P1k(state)+2100;
      
        H2_a(state) = refpropm('H','P',Pchoke_a(state),'S',S_a(state),'CO2');
        V2_a(state) = sqrt(2*(H_a(state)-H2_a(state)));
    end
    
    % Guess initial choke pressure
    Pchoke_b(state) = min(13000, P2k(state)-1000) ;% kPa
    Ma_b(state) = .5;
    while abs(Ma_b(state)-1) > .001
        SS_b(state) = refpropm('A','P',Pchoke_b(state),'S',S_b(state),'CO2');
        H2_b(state) = refpropm('H','P',Pchoke_b(state),'S',S_b(state),'CO2');
        V2_b(state) = sqrt(2*(H_b(state)-H2_b(state)));
        Ma_b(state) = V2_b(state)/SS_b(state);
        Pchoke_b(state) = Pchoke_b(state) + (Ma_b(state) - 1) * 50;
    end
    % Final pressure = P1 + delta P
    P3(state) = P1k(state) + 2.1; %MPa
    
    H3_a(state) = H_a(state);
    H3_b(state) = H_b(state);
    H_final(state) = refpropm('H','T',150+273.15,'P',P3(state),'CO2');
     
    rho2_a(state) = refpropm('D', 'P', Pchoke_a(state), 'S', S_a(state), 'CO2');
    rho2_b(state) = refpropm('D', 'P', Pchoke_b(state), 'S', S_b(state), 'CO2');
    rho3_a = refpropm('D','H', H3_a(state), 'P', P3(state),'CO2') ;
    
    T3_a =  refpropm('T','H', H3_a(state), 'P', P3(state), 'CO2');
    
    % Calculate mass flow by rearranging the following equations
    %Md_A + md_B = .26
    %Md_A * H_a + Md_B * H_b = H_final(state) * .26
    Md_B(state) = .26*(H_final(state)-H_a(state))/(H3_b(state)-H3_a(state));
    Md_A(state) = .26 - Md_B(state);
    
    % Calculate Valve area using mdot = rho*A*V
    Area_a(state)= Md_A(state)./(V2_a(state)*rho2_a(state))
    Area_b(state)= Md_B(state)./(V2_b(state)*rho2_b(state))
    
    D_a(state)= sqrt(4*Area_a(state)/pi);
    D_b(state)= sqrt(4*Area_b(state)/pi);

    Cd=.61;
    Cv_a(state)= 46250.9*Cd*D_a(state)^2
    Cv_b(state)= 46250.9*Cd*D_b(state)^2
end
%%
% Efficiency loss from siphoning mass from the system
ef_m = (m_dot - .26)./m_dot 

% Finds maximum power based on actual mass flow rate in the state
p_actual = 400*(m_dot/6.5)

% Accounting for power loss due to valves
p_actuator = (24 * 20/1000)/1000

% Accounting for power loss due to computer for valve control 
p_computer = .1 % average power of a computer, kWatts 

% Overall efficiency with valves and computer accounted for
eff = ef_m .* (p_actual-(p_actuator + p_computer))./p_actual

% Power output based on efficiency 
p_output = eff .* p_actual 

%%
figure
title('Enthalpy and Mass Flow Rate')
yyaxis left
y = [Md_A(1:9); Md_B(1:9)]'

ba1 = bar(y, .5 ,'stacked','FaceColor','flat');
ylabel('Mass Flow Rate (Kg/s)', 'Color' ,'k')
ba1(1).CData = [ 0 .4 1];
ba1(2).CData = [1 .4 0];

yyaxis right
plot(1:9, H_a(1:9), 'Color', [0/256 0/256 245/256], 'MarkerSize', 4, 'LineWidth',3)
hold on

plot(1:9, H_b(1:9), 'Color', [1,0,0], 'MarkerSize', 4, 'LineWidth',3)
hold on 
p = plot(1:9, H_final(1:9),'Color', [27/256 1/256 111/256],'MarkerSize' , 4, 'LineWidth',3)
legend('Mdot 3', 'Mdot 4','H3', 'H4', 'H final')

xticklabels({'1: T=-20 ', '2: T=-10 ', '3: T=0 ', '4: T=10 ', '5: T=20 ', '6: T=30 ' ...
    '7: T=40 ', '8: T=50 ', 'Startup: T=20 '})
xtickangle(45)
xlabel('State (Ambient Condition)')
ylabel('Enthalpy (kJ/kG)',  'Color' ,'k')
