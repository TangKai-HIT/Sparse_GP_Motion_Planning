%demo1: construct sparse 1-D GP trajectory
close all; clear; clc;

%% Init min-jerk sparse GP trajectory
% construct anchor states using cubic-spline
t_wp = [0, 5, 7, 8, 10, 15, 18];
q_wp = [3; -2; -5; 0; 6; 12; 8];
init_vel = 0;
end_vel = 0;
numSample = 10;

[q, dq, ddq, dddq, t_sample, ~]=cubicSplineTraj_wp(q_wp, t_wp, init_vel, end_vel, numSample);

%origin spline plot points
[q_origin, ~, ~, ~, t_origin, ~]=cubicSplineTraj_wp(q_wp, t_wp, init_vel, end_vel, 500);

% init prior params
Q_c = 1e-2; %noise variance at each dimension
kappa0 = diag([1e-1, 1e-1, 1e-1]); %init state variance
supportPts_minjerk = [q, dq, ddq];

minJerkTraj = GPTrajSparseMinJerk(supportPts_minjerk, t_sample, Q_c, supportPts_minjerk(1, :)', kappa0);

minJerkTraj.updateLifted();

%% Plot min-jerk sparse GP trajectory
figure('Position',[500,100, 1000,500]);
plotSparseGP_1D(gca, minJerkTraj);
title("sparse GP min-jerk prior interpolation");
hold on;
plot(t_origin, q_origin, '--b', 'LineWidth', 0.8);

%% Init min-acc sparse GP trajectory
% init prior params
Q_c = 1e-1; %noise variance at each dimension
kappa0 = 1e-1 * eye(2); %init state variance
supportPts_minacc = [q, dq];

minAccTraj = GPTrajSparseMinAcc(supportPts_minacc, t_sample, Q_c, supportPts_minacc(1, :)', kappa0);

minAccTraj.updateLifted();

%% Plot min-acc sparse GP trajectory
figure('Position',[500,100, 1000,500]);
plotSparseGP_1D(gca, minAccTraj);
title("sparse GP min-acc prior interpolation");
hold on;
plot(t_origin, q_origin, '--b', 'LineWidth', 0.8);