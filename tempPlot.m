
p2 = 26e6
md = 9
k = 1e-12 ;% 
A = 20 ;% m^-1
L = 112 ;% m
Pa = 32.5e6 ;
Tb = 129 + 273; % Temp at bottom

[p23, friction_loss_down, gravity_gain_down, temp_down] = pressure_drop_down(md, p2, 303);
p3 = p23(end)
density = refpropm('D','T',Tb,'P',p3/ 1e3, 'CO2');
dynamic_v = refpropm('V','T',Tb,'P',p3/1e3, 'CO2');

t3 = temp_down(end);

final_mdot = md - 3.24;
Q = final_mdot/ density;
delta_p = Q * dynamic_v* L/(k*A);
t4 = 129+273;
p4 = p3 - delta_p;
[p5, friction_loss_up, gravity_loss_up, temp_top] = pressure_drop_up2(final_mdot, p4);   
t5 = temp_top;


temp = [];
temp1 = [];
temp2 = [];
x = [];
h = [];
x1 = [];
h1 = [];
x2 = [];
h2 = [];
delta_l = 1;
i =  1
for height = 3200 - delta_l : -1 * delta_l : 0
   temp = [temp_top(i) temp];
   x = [112 x];
   h = [height h ];
   i = i+1;
end
i = length(temp_down)
for height = 0 : delta_l : 3200 - delta_l      
   temp1 = [temp_down(i) temp1];
   h1 = [h1 height];
   x1 =[0 x1];
   i = i-1
end

for i = 1:112
    x2 = [i x2];
    h2 = [3400 h2];
    temp2 = [129 + 273 temp2];

end

%%
figure

ax1 = axes;
hh1 = scatter(ax1, x1 ,h1 , [],temp1, 'filled' , 's')
ax1.YDir = 'reverse';
hold on


hh2 = scatter(ax1 ,x2 ,h2 , [], 'red', 'filled' , 's')

%yyaxis right
ax2 = axes;

hh = scatter(ax2, x ,h , [],temp, 'filled' , 's')
ax2.YDir = 'reverse'
title('Temperature Gradient of CO2')
xlabel('Distance Traveled in X (m)')
ylabel('Depth (m)')
%axes([0 3200 min(temp) max(temp)])
%ah1 = axes;
%%Link them together
linkaxes([ax1,ax2])
%%Hide the top axes
ax2.Visible = 'off';

ax1.Visible = 'on';
ax2.XTick = [];
ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'copper')
colormap(ax2,'autumn')
grid on
%%Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.07 .11 .03 .815]);
cb2 = colorbar(ax2,'Position',[.90 .11 .03 .815]);
title('Temperature Gradient of CO2')
xlabel('Distance Traveled in X (m)')
ylabel('Depth (m)')
grid on
s = 25; %Marker width in units of X
%Obtain the axes size (in axpos) in Points
currentunits = get(gca,'Units');
set(gca, 'Units', 'Points');
axpos = get(gca,'Position');
set(gca, 'Units', currentunits);
markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
set(hh, 'SizeData', markerWidth^2)
set(hh1, 'SizeData', markerWidth^2)
title('Temperature Gradient of CO2')
xlabel('Distance Traveled in X (m)')
ylabel('Depth (m)')
%set(hh2, 'SizeData', markerWidth^2)
grid on