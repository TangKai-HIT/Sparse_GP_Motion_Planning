%%GPMP 2D point: use gpmp planner to plan a 2D point trajectory with fixed boundary points
clc; clear; close all;

%% Construct 2D cost map
bbox = [0 1 0 1]; %unit bounding box
% rectangle 1
rectangle_array(1).low = [0.2 0.2];
rectangle_array(1).high = [0.4 0.5];
% rectangle 2
rectangle_array(2).low = [0.6 0.5];
rectangle_array(2).high = [0.8 0.9];
% Get map
resolution = 0.001; %resolution of map
map = rectangle_maps( bbox, rectangle_array, resolution);
% Cost map 
epsilon = 0.12; % horizon upto which distance is computed
cost_map = create_costmap_sqdist(map, epsilon);
% Cost map derivatives
[cost_map_x, cost_map_y] = get_cost_map_derivatives(cost_map);
% Visualize cost map & robot
figure;
visualize_cost_map(cost_map);
visualize_rectangles(gca, rectangle_array, 'r', '--');

%% Define & init 2D point support states
p_start = [0 0.5];
p_end = [1 0.5];

% cubic traj samples (cubic traj is closer to the min-jerk prior(3-order))
N = 10;
t_support = linspace(0, 3, N);
[p_sup, pd_sup, pdd_sup]=cubicpolytraj([p_start', p_end'], t_support([1, N]), t_support);
supportPts = [p_sup', pd_sup', pdd_sup'];

% trapveltraj samples 
% (note: don't use non smooth curve as initiation, otherwise the states would move far away from original place)
% N = 5;
% [p_sup, pd_sup, pdd_sup, t_support, ~]=trapveltraj([p_start', p_end'], N,"PeakVelocity",0.5);
% % hold on; plot(p(1, :), p(2, :), '-b', 'LineWidth', 1); 
% supportPts = [p_sup', pd_sup', pdd_sup'];

%% Define params, cost functions & init gpmp planner
% init prior params
Q_c = diag([0.1, 0.1]); %noise variance at each dimension
kappa0 = diag(1e-1*ones(1, 2*3)); %init state variance

% init coeffs
omega = 3000; %obs cost coeffs
lambda = 1; %prior cost coeffs
eta = 100; %regularization coeffs

% set options
options = gpmp_optimset_init();
options.FixedStateId = [1, 5];
options.MaxIter = 500;
options.TolFun = 1e-2;

% init sparse gp traj & gpmp planner
gpSparseSets = gpSparse_init_set(supportPts, t_support, Q_c, supportPts(1, :)', kappa0);
% minJerkTraj = GPTrajSparseMinJerk(supportPts, t_support, Q_c, supportPts(1, :)', kappa0);
varDim = 2; %2D

gpmpPlanner2D = gpmpBound2D(gpSparseSets, varDim, eta, omega, lambda, options);
% gpmpPlanner2D = gpmpUncon2D(minJerkTraj, varDim, order, eta, lambda, options);

% SDF function handles
gpmpPlanner2D.disFieldfn = @(xi) value_wpset_map( xi, cost_map ); % distance field evaluation function handle
gpmpPlanner2D.disGradfn = @(xi) [value_wpset_map( xi, cost_map_x ), value_wpset_map( xi, cost_map_y )]; % distance field gradient evaluation function handle

% upsampling
resolution = 0.05;
gpmpPlanner2D.upSampleByResol(resolution * ones(1, gpmpPlanner2D.numIntervals));

%% Perform gpmp
tic;
results = gpmpPlanner2D.solve();
solved_time = toc; 

fprintf("solve time: %.4f (s)\n", solved_time);
fprintf("iterations: %d\n", results.steps);
fprintf("exitFlag: %d\n", results.exitFlag);

%% Visualize results
gpmpPlanner2D.plotDeformHis(gca, @winter , 10);
gpmpPlanner2D.plotCostHistory();

%% Plot velcity, acceleration, jerk
%plot
gpmpPlanner2D.gpTrajSparse.supportPts = results.x_solve;
[plotStates, plotId] = sampleStatesHelper(gpmpPlanner2D.gpTrajSparse, 500);
diff1Mat = getDiffMatrix_quadInterp(plotId, 1);

figure("Position",[500,100, 800, 800]);

subplot(3,1,1);
plot(plotId, plotStates(:, [3,4])); legend('x', 'y');
title("vel")
subplot(3,1,2);
plot(plotId, plotStates(:, [5,6])); legend('x', 'y');
title("accel")
subplot(3,1,3);
plot(plotId, diff1Mat*plotStates(:, [5,6])); legend('x', 'y');
title("jerk")