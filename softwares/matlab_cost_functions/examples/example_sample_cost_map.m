%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%% EXAMPLE_SAMPLE_COST_MAP
% Create sample cost map based on map created on rectangles
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

%% Visualize map
figure;
visualize_map(map);

%% Cost map 
epsilon = 0.15; % horizon of expansion
cost_map = create_costmap_sqdist(map, epsilon);

%% Visualize cost map
figure;
visualize_cost_map(cost_map);

%% Cost map derivatives
[cmap_x, cmap_y] = get_cost_map_derivatives(cost_map);
figure;
visualize_cost_map(cmap_x);
figure;
visualize_cost_map(cmap_y);
