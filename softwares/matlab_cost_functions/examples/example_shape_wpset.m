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

% Visualize 
figure;
visualize_shapes(rectangle_array);
axis(bbox)

%% Create a trajectory
p_start = [0.25 0.7];
p_end = [0.7 0.25];
N = 20;
wpset = get_wpset(p_start, p_end, N);

% snip of end points
wpset = wpset(2:(end-1),:);

%% Check cost and gradient
epsilon = 0.15;
c_obs_fn_wp = @(pt) costfn_sqdist_value_pt_shapes(pt, rectangle_array, epsilon);
grad_c_obs_fn_wp = @(pt) costfn_sqdist_grad_pt_shapes(pt, rectangle_array, epsilon);

c_obs_fn_arclength_wpset = @(xi) value_arclength_fn_wp( xi, c_obs_fn_wp, p_start );
grad_c_obs_fn_arclength_wpset = @(xi) gradient_arclength_fn_wp( xi, c_obs_fn_wp, grad_c_obs_fn_wp, p_start);

cost = c_obs_fn_arclength_wpset(wpset)
grad = grad_c_obs_fn_arclength_wpset(wpset)
