%%GPMP 2D states (scene 2): use gpmp planner to plan a 2D state trajectory with floating boundary states (with controllable variance)
clc; clear; close all;

%% Construct 2D cost map
bbox = [0 1 0 1]; %unit bounding box
% rectangle 1
rectangle_array(1).low = [0.25 0.3];
rectangle_array(1).high = [0.45 0.5];
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
p_start = [0.1 0.7];
p_end = [1 0.5];

dp_start = [0 -0.2];
dp_end = [0 0];

% cubic traj samples (cubic traj is closer to the min-jerk prior(3-order))
N = 10;
t_support = linspace(0, 3, N);
[p_sup, pd_sup, pdd_sup]=cubicpolytraj([p_start', p_end'], t_support([1, N]), t_support, "VelocityBoundaryCondition",[dp_start; dp_end]');
supportPts = [p_sup', pd_sup', pdd_sup'];

%% Define params, cost functions & init gpmp planner
% init prior params
Q_c = diag([0.1, 0.1]); %noise variance at each dimension
kappa0 = diag(0.5*1e-3*ones(1, 2*3)); %init state variance
kappa_goal = kappa0;

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
order = 3;

gpmpPlanner2D = gpmpFloatBound2D(gpSparseSets, varDim, order, eta, omega, lambda, kappa_goal, options, []);

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

%% Plot variance
t_samp = linspace(t_support(1), t_support(end), 100);
[~, ~, ~, ~, observedPts_all, kappa_all, ~]  = gpmpPlanner2D.gpTrajSparse.upSampleInterval_observeGoal(1, N, t_samp(2:end-1), kappa0, kappa_goal);

kappa_main = diag(kappa_all);
x_Id = 1 : 6 : 100*6;
y_Id = 2 : 6 : 100*6;
varScale = 0.5e1;
patch(gca, [observedPts_all(x_Id); observedPts_all(flip(x_Id))], ...
		[observedPts_all(y_Id)+squeeze(abs(kappa_main(y_Id)).^.5)*varScale; observedPts_all(flip(y_Id))-squeeze(abs(kappa_main(flip(y_Id))).^.5)*varScale], ...
		[.6 .8 .8],'edgecolor','none'); %added abs to avoid negative value around 0 (numerical problem) 
alpha 0.6

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