%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%% EXAMPLE_COSTFN_GEOMETRIC_LINE_MAP
% Interactive example to check geometric line cost in a map
clc;
clear;
close all;

%% Create two rectangles
bbox = [0 1 0 1]; %unit bounding box
% rectangle 1
rectangle_array(1).low = [0.2 0.2];
rectangle_array(1).high = [0.4 0.6];
% rectangle 2
rectangle_array(2).low = [0.6 0.5];
rectangle_array(2).high = [0.8 0.9];
resolution = 0.001; %resolution of map

% Get map
map = rectangle_maps( bbox, rectangle_array, resolution);


%% Visualize 
figure;
visualize_map(map);
grid on;
hold on;

%% Interactive session to draw line and see cost
infinity = inf; %cost of a line which is in collision
while(1)
    title('Click two points to create a line');
    [x, y] = ginput(2);
    plot(x, y, 'b');
    c = costfn_geometric_value_line_map([x(1) y(1)], [x(2) y(2)], map, infinity);
    str = sprintf('Cost is %f , Hit return to continue', c);
    title(str);
    pause;
    clf;
    visualize_map(map);
    grid on;
    hold on;
end
