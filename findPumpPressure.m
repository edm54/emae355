% Finds the pressure needed at the top of pipe to reach a given
% pressure needed at the bottom

p4 = 0;

% Initial Guess
p3 = 15e6;

% Mass flow rate 
md = 3.24;
md = 10
while (32.5e6 - p4 > 1e3)
    [p4, press_drop_down, gravity_gain_down, temp_down] = pressure_drop_down(md, p3);
    p3 = p3 + .5 * (32.5e6 - p4);
end

disp('Pressure needed:')
disp(p3)
