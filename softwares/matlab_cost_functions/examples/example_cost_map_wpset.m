%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%% EXAMPLE_COST_MAP_WPSET
% Create sample cost map and then do operations on wpset
clc;
clear;
close all;

%% Create cost map
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

% Cost map 
epsilon = 0.15; % horizon of expansion
cost_map = create_costmap_sqdist(map, epsilon);

% Cost map derivatives
[cost_map_x, cost_map_y] = get_cost_map_derivatives(cost_map);

% Visualize cost map
figure;
visualize_cost_map(cost_map);

%% Create a trajectory
p_start = [0.25 0.7];
p_end = [0.7 0.25];
N = 20;
wpset = get_wpset(p_start, p_end, N);

% snip of end points
wpset = wpset(2:(end-1),:);

%% Check cost and gradient
cost = value_arclength_wpset_map( wpset, cost_map, p_start )
grad = gradient_arclength_wpset_map( wpset, cost_map, cost_map_x, cost_map_y, p_start)

% check using another way
c_obs_fn_wp = @(pt) value_wpset_map(pt, cost_map);
grad_c_obs_fn_wp = @(pt) [value_wpset_map(pt, cost_map_x) value_wpset_map(pt, cost_map_y)];

c_obs_fn_arclength_wpset = @(xi) value_arclength_fn_wp( xi, c_obs_fn_wp, p_start );
grad_c_obs_fn_arclength_wpset = @(xi) gradient_arclength_fn_wp( xi, c_obs_fn_wp, grad_c_obs_fn_wp, p_start);

cost2 = c_obs_fn_arclength_wpset(wpset)
grad2 = grad_c_obs_fn_arclength_wpset(wpset)
