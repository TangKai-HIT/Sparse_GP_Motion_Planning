%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%% EXAMPLE_COSTFN_GEOMETRIC_LINE_SHAPE
% Interactive example to check geometric line cost
clc;
clear;
close all;

%% Create two rectangles
bbox = [0 1 0 1]; %unit bounding box
% rectangle 1
rectangle_array(1) = get_rectangle_shape(0.2,0.2,0.2,0.4);
% rectangle 2
rectangle_array(2) = get_rectangle_shape(0.6,0.5,0.2,0.4);

%% Visualize 
figure;
visualize_shapes(rectangle_array);
axis(bbox)
grid on;

%% Interactive session to draw line and see cost
infinity = inf; %cost of a line which is in collision
while(1)
    title('Click two points to create a line');
    [x, y] = ginput(2);
    plot(x, y, 'b');
    c = costfn_geometric_value_line_shapes([x(1) y(1)], [x(2) y(2)], rectangle_array, infinity);
    str = sprintf('Cost is %f , Hit return to continue', c);
    title(str);
    pause;
    clf;
    hold on;
    visualize_shapes(rectangle_array);
    axis(bbox);
    grid on;
end
