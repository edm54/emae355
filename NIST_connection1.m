%%
%load(refpropm)
% First Value is startup condition
%9 is a copy of 5 
Ambient= [ -20, -10, 0, 10, 20, 30, 40, 50, 20, 20];
T1 = [ -20,-10, 0 ,10, 20, 30, 40, 50, 20, 20];
T2 = [-14, -2, 10,23, 38, 57, 71, 85, 27, 140];
T3 = [ 33,39,47, 59, 71,85,104, 122, 79, 90];
T4 = [ 200,200,200,200,200,200,200,200, 200, 210];
T5 = [81, 81,87, 96, 103, 113, 133, 150, 184, 103];
T6 = [ -2, 3, 15, 28,43, 62, 77, 91, 33, 43];
P1 = [ 3.3,3.3, 3.9, 5, 6.2, 7.4, 9.2, 10, 6.2, 6.2];
P2 = [ 14, 15, 16, 18, 20, 20.7, 19.9, 17.8, 11, 20] ;
m_dot = [ 4.9, 4.9, 5.2, 5.8, 6.4, 6.5, 5.5, 4, 3, 6.4];

T1K = T1+ 273.15;
T2K = T2 + 273.15;
T3K = T3+ 273.15;
T4K = T4 + 273.15;
T5K = T5+ 273.15;
T6K = T6 + 273.15;
P1k = P1 * 1000;
P2k = P2 * 1000;


%pressure = p2
%H = refpropm('H','T',373.15,'P', 1,'CO2');
%specfic_v = 1/refpropm('D','T', 373.15, 'P',1, 'CO2');
%S = refpropm('S','T',373.15,'P',1,'CO2');

for state = 1:length(T1)
    H_a(state) = refpropm('H','T',T3K(state),'P',P2k(state),'CO2');
    specfic_v_a(state) = 1/refpropm('D','T',T3K(state),'P',P2k(state),'CO2');
    S_a(state) = refpropm('S','T',T3K(state),'P',P2k(state),'CO2');
    
    H_b(state) = refpropm('H','T',T4K(state),'P',P2k(state),'CO2');
    specfic_v_b(state) = 1/refpropm('D','T',T4K(state),'P',P2k(state),'CO2');
    S_b(state) = refpropm('S','T',T4K(state),'P',P2k(state),'CO2');  
end
% u = mdot/rho* a 
% can input a 
% can get c
% need m dot until u = c
for state = 1:length(T1)
%or state = 
    try
    Pchoke_a(state) = min(13000, P2k(state)-1000);% kPa
    Ma_a(state) = .5;
        while abs(Ma_a(state)-1) > .001
            %disp("p "+ Pchoke_a(state));
            %disp("s "+ S_a(state))
            SS_a(state) = refpropm('A','P',Pchoke_a(state),'S',S_a(state),'CO2');
            %disp("ss "+ SS_a(state));
            H2_a(state) = refpropm('H','P',Pchoke_a(state),'S',S_a(state),'CO2');
            %disp("H2a "+ H2_a(state));

            V2_a(state) = sqrt(2*(H_a(state)-H2_a(state)));
            %disp("V2 " + V2_a(state));
            Ma_a(state) = V2_a(state)/SS_a(state);
            Pchoke_a(state) = Pchoke_a(state) + (Ma_a(state) - 1) * 50;
            %disp("P " + Pchoke_a(state))
            %disp("Ma " +   Ma_a(state));
        end
    catch ME
        disp(ME)
        %Sets the choke pressure equal to the final pressure when the mach
        %state can not be reached. This is because the flow pressure can
        %now communicate up to the choke without the Mach 1 flow
        Pchoke_a(state) = P1k(state)+2100;
        
        %Same as above
        %SS_a(state) = refpropm('A','P',Pchoke_a(state),'S',S_a(state),'CO2');
        %disp("ss "+ SS_a(state));
        H2_a(state) = refpropm('H','P',Pchoke_a(state),'S',S_a(state),'CO2');
        %disp("H2a "+ H2_a(state));

        V2_a(state) = sqrt(2*(H_a(state)-H2_a(state)));
        %disp("V2 " + V2_a(state));
    end
    %disp("a")
    %disp(V2_a(state))
    %disp(H2_a(state))
    %disp(SS_a(state))
    %disp(Ma_a(state))
    Pchoke_b(state) = min(13000, P2k(state)-1000) ;% kPa
    Ma_b(state) = .5;
    while abs(Ma_b(state)-1) > .001
        %disp("p "+ Pchoke_b(state))
        %disp("s "+ S_b(state))
        SS_b(state) = refpropm('A','P',Pchoke_b(state),'S',S_b(state),'CO2');
        %disp("ss "+ SS_b(state));
        H2_b(state) = refpropm('H','P',Pchoke_b(state),'S',S_b(state),'CO2');
        %disp("H2 "+ H2_b(state));
        V2_b(state) = sqrt(2*(H_b(state)-H2_b(state)));
        %disp("V2 " + V2_b(state));
        Ma_b(state) = V2_b(state)/SS_b(state);
        Pchoke_b(state) = Pchoke_b(state) + (Ma_b(state) - 1) * 50;
    end
    
    %disp(V2_b(state));
    %disp(H2_b(state));
    %disp(SS_b(state));
        
    rho2_a(state) = refpropm('D', 'P', Pchoke_a(state), 'S', S_a(state), 'CO2');
    rho2_b(state) = refpropm('D', 'P', Pchoke_b(state), 'S', S_b(state), 'CO2');
    
    H3_a(state) = H_a(state);
    H3_b(state) = H_b(state);
    P3(state) = P1k(state) + 2.1; %MPa
   
    rho3_a = refpropm('D','H', H3_a(state), 'P', P3(state),'CO2') ;   
    T3_a =  refpropm('T','H', H3_a(state), 'P', P3(state), 'CO2');
    
    H_final(state) = refpropm('H','T',150+273.15,'P',P3(state),'CO2');
    
    %energy equation
    %so vary a, which will vry m dot which is constrainted by 1+2= .26 
    % m dot is related to thre choke area
    area = 0;
    Md_B(state) = .26*(H_final(state)-H_a(state))/(H3_b(state)-H3_a(state))
    Md_A(state) = .26 - Md_B(state)
    
    %Md_A + md_B = .26
    %Md_A * H_a + Md_B * H_b = H_final(state) * .26
    
    %get md_a and md_b
    %mdot = rho v a
    %v at choke, rho at chock, get a 
    
    % T = 150 C = 423.15 K 
    %while 
    disp(state)
    
    Area_a(state)= Md_A(state)./(V2_a(state)*rho2_a(state))
    Area_b(state)= Md_B(state)./(V2_b(state)*rho2_b(state))
    
    D_a(state)= sqrt(4*Area_a(state)/pi);
    D_b(state)= sqrt(4*Area_b(state)/pi);

    Cd=.61;
    Cv_a(state)= 46250.9*Cd*D_a(state)^2
    Cv_b(state)= 46250.9*Cd*D_b(state)^2
end
%%
figure
title('Enthalpy and Mass Flow Rate')

yyaxis left
y = [Md_A(1:8); Md_B(1:8)]'

ba = bar(y, .5 ,'stacked','FaceColor','flat')
ylabel('Mass Flow Rate (Kg/s)', 'Color' ,'k')
ba(1).CData = [ 0 .4 1];
ba(2).CData = [1 .4 0];

yyaxis right
plot(1:8, H_a(1:8), 'Color', [14/256 5/256 51/256], 'MarkerSize', 4, 'LineWidth',3)
hold on 
plot(1:8, H_b(1:8), 'Color', [56/256 9/256 209/256], 'MarkerSize', 4, 'LineWidth',3)
hold on 
p = plot(1:8, H_final(1:8),'Color', [27/256 1/256 111/256],'MarkerSize' , 4, 'LineWidth',3)
legend('Mdot A', 'Mdot B','H3', 'H4', 'H final')

%set(p, {'color'}, {[115/256 75/256 247/256]; [56/256 9/256 209/256] ; [27/256 1/256 111/256]});
%colors = hsv(5);
%set(p, {'color'}, num2cell(colors, 5));
xlabel('State')
ylabel('Enthalpy (kJ/kG)',  'Color' ,'k')


%yylabel('Enthalpy')