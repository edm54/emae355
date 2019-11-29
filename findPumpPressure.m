% Finds the pressure needed at the top of pipe to reach a given
% pressure needed at the bottom
p3 = 0;

% Initial Guess
p2 = 15e6;

% Mass flow rate 
md = 3.24;
%md = 10
while (32.5e6 - p3 > 1e3)
    [p3, press_drop_down, gravity_gain_down, temp_down] = pressure_drop_down(md, p2, 30+273);
    p2 = p2 + .5 * (32.5e6 - p3);
end

disp('Pressure needed:')
disp(p2)
