%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%% EXAMPLE_COSTFN_SMOOTH
% Example showing smoothness cost function
clc;
clear;
close all;

%% Create a trajectory
radius = 1;
N = 20;
thet = transpose(linspace(0, pi, N));
wpset = [radius*cos(thet) radius*sin(thet)];

p_start = wpset(1,:);
p_end = wpset(end,:);
wpset = wpset(2:(end-1),:);

%% Smoothness cost
[ A, b, c, K, e ] = create_smooth_matrices( p_start, p_end, size(wpset,1) );
cost = costfn_smooth_value_wpset(wpset, A, b, c)
grad = costfn_smooth_grad_wpset(wpset, A, b)

%% Plot
figure;
hold on;
xi = wpset;
for i = 1:10
    cost = costfn_smooth_value_wpset(xi, A, b, c);
    grad = costfn_smooth_grad_wpset(xi, A, b);
    fprintf('Iter: %d Cost: %f\n', i, cost);
    xi_disp = [p_start; xi; p_end];
    plot(xi_disp(:,1), xi_disp(:,2));
    xi = xi - 0.01*(A\grad);
end
