%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%% EXAMPLE_POINT_SQDIST_SHAPE
% Interactive example to check sqdist cost of a point
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
epsilon = 0.15;
while(1)
    title('Click a point to check its sqdist cost and grad');
    [x, y] = ginput(1);
    [c, c_grad] = costfn_sqdist_pt_shapes([x y], rectangle_array, epsilon);
    str = sprintf('Cost is %f Grad is %f, %f, Hit return to continue', c, c_grad(1), c_grad(2));
    title(str);
    pause;
    clf;
    hold on;
    visualize_shapes(rectangle_array);
    axis(bbox);
    grid on;
end
