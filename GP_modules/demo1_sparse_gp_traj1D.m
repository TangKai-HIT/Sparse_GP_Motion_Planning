%demo1: construct sparse 1-D GP trajectory
close all; clear; clc;

%% Init min-jerk sparse GP trajectory
% construct anchor states using cubic-spline
t_wp = [0, 5, 7, 8, 10, 15, 18];
q_wp = [3; -2; -5; 0; 6; 12; 8];
init_vel = 0;
end_vel = 0;
numSample = 10;

[q_cubic, dq_cubic, ddq_cubic, dddq_cubic, t_sample, ~]=cubicSplineTraj_wp(q_wp, t_wp, init_vel, end_vel, numSample);

%origin spline plot points
[q_origin, dq_origin, ddq_origin, dddq_origin, t_origin, ~]=cubicSplineTraj_wp(q_wp, t_wp, init_vel, end_vel, 500);

% init prior params
Q_c = 1e-2; %noise variance at each dimension
kappa0 = diag([1e-1, 1e-1, 1e-1]); %init state variance
supportPts_minjerk = [q_cubic, dq_cubic, ddq_cubic];

minJerkTraj = GPTrajSparseMinJerk(supportPts_minjerk, t_sample, Q_c, supportPts_minjerk(1, :)', kappa0);

minJerkTraj.updateLifted();

%% Plot min-jerk sparse GP trajectory
figure('Position',[500,80, 1000,900]);
subplot(3,1,1);
[minJerkTrajStates, sampleId_minJerk] = plotSparseGP_1D(gca, minJerkTraj, 1, 500);
title("position interpolation");
hold on;
plot(t_origin, q_origin, '--b', 'LineWidth', 0.8);
legend("conditioned variance","mean prediction", "support states", "ground truth", "Location","northwest");

subplot(3,1,2);
plotSparseGP_1D(gca, minJerkTraj, 2, 500);
title("velocity interpolation");
hold on;
plot(t_origin, dq_origin, '--b', 'LineWidth', 0.8);
legend("conditioned variance","mean prediction", "support states", "ground truth", "Location","northwest");

subplot(3,1,3);
plotSparseGP_1D(gca, minJerkTraj, 3, 500);
title("acceleration interpolation");
hold on;
plot(t_origin, ddq_origin, '--b', 'LineWidth', 0.8);
legend("conditioned variance","mean prediction", "support states", "ground truth", "Location","northwest");

sgtitle("Sparse GP min-jerk prior interpolation");

%% Init min-acc sparse GP trajectory
% init prior params
Q_c = 1e-1; %noise variance at each dimension
kappa0 = 1e-1 * eye(2); %init state variance
supportPts_minacc = [q_cubic, dq_cubic];

minAccTraj = GPTrajSparseMinAcc(supportPts_minacc, t_sample, Q_c, supportPts_minacc(1, :)', kappa0);

minAccTraj.updateLifted();

%% Plot min-acc sparse GP trajectory
figure('Position',[500,80, 1000,800]);
subplot(2,1,1)
[minAccTrajStates, sampleId_minAcc] = plotSparseGP_1D(gca, minAccTraj, 1, 500);
title("position interpolation");
hold on;
plot(t_origin, q_origin, '--b', 'LineWidth', 0.8);
legend("conditioned variance","mean prediction", "support states", "ground truth", "Location","northwest");

subplot(2,1,2)
plotSparseGP_1D(gca, minAccTraj, 2, 500);
title("velocity interpolation");
hold on;
plot(t_origin, dq_origin, '--b', 'LineWidth', 0.8);
legend("conditioned variance","mean prediction", "support states", "ground truth", "Location","northwest");

sgtitle("Sparse GP min-acc prior interpolation");

%% Test acceleration cost
cubic_spline_cost = dot(sum(ddq_origin(1:end-1) .* ddq_origin(1:end-1), 2), diff(t_origin));
disp("cubic spline acceleration cost: "); disp(cubic_spline_cost);

minAccTraj_acc = diff(minAccTrajStates(:, 2)) ./ diff(sampleId_minAcc);
minAccTraj_cost = dot(sum(minAccTraj_acc .* minAccTraj_acc, 2), diff(sampleId_minAcc));
disp("min-acc sparse GP acceleration cost: "); disp(minAccTraj_cost);

minSnapTraj_jerkCost = dot(sum(minJerkTrajStates(1:end-1, 3) .* minJerkTrajStates(1:end-1, 3), 2), diff(sampleId_minJerk));
disp("min-jerk sparse GP acceleration cost: "); disp(minSnapTraj_jerkCost);

%% Init min-snap sparse GP trajectory
t_wp = [0,  10, 20,  30];
q_wp = [-5;  0;  8;  20;];

t_Samples = linspace(t_wp(1), t_wp(end), 500);
[q_quin, dq_quin, ddq_quin, ~] = quinticpolytraj(q_wp', t_wp, t_Samples);
dddq_quin = diff(ddq_quin) ./ diff(t_Samples);

id_sq = [1: 50: 450, 500];
q_quin_sup = q_quin(id_sq);
dq_quin_sup = dq_quin(id_sq);
ddq_quin_sup = ddq_quin(id_sq);
dddq_quin_sup = dddq_quin([1: 50: 450, 500-1]);
t_Support = t_Samples(id_sq);

% init prior params
Q_c = 1e-2; %noise variance at each dimension
kappa0 = diag([1e-1, 1e-1, 1e-2, 1e-2]); %init state variance
supportPts_minSnap = [q_quin_sup; dq_quin_sup; ddq_quin_sup; dddq_quin_sup]';

minSnapTraj = GPTrajSparseMinSnap(supportPts_minSnap, t_Support, Q_c, supportPts_minSnap(1, :)', kappa0);

minSnapTraj.updateLifted();

%% Plot min-snap sparse GP trajectory
figure('Position',[100,80, 1500,900]);
subplot(2,2,1);
[minSnapTrajStates, sampleId_minSnap] = plotSparseGP_1D(gca, minSnapTraj, 1, 500 , 5e1);
title("position interpolation");
hold on;
plot(t_Samples, q_quin, '--b', 'LineWidth', 0.8);
legend("conditioned variance","mean prediction", "support states", "ground truth", "Location","northwest");

subplot(2,2,2);
plotSparseGP_1D(gca, minSnapTraj, 2, 500, 2.5e1);
title("velocity interpolation");
hold on;
plot(t_Samples, dq_quin, '--b', 'LineWidth', 0.8);
legend("conditioned variance","mean prediction", "support states", "ground truth", "Location","northwest");

subplot(2,2,3);
plotSparseGP_1D(gca, minSnapTraj, 3, 500, 2.5e1);
title("acceleration interpolation");
hold on;
plot(t_Samples, ddq_quin, '--b', 'LineWidth', 0.8);
legend("conditioned variance","mean prediction", "support states", "ground truth", "Location","northwest");

subplot(2,2,4);
plotSparseGP_1D(gca, minSnapTraj, 4, 500, 2.5e1);
title("jerk interpolation");
hold on;
plot(t_Samples(1:end-1), dddq_quin, '--b', 'LineWidth', 0.8);
legend("conditioned variance","mean prediction", "support states", "ground truth", "Location","northwest");

sgtitle("Sparse GP min-Snap prior interpolation");

%% Test jerk cost
quintic_jerkCost = dot(dddq_quin .* dddq_quin, diff(t_Samples));
disp("quintic trajectory jerk cost: "); disp(quintic_jerkCost);

minSnapTraj_jerkCost = dot(sum(minSnapTrajStates(1:end-1, 4) .* minSnapTrajStates(1:end-1, 4), 2), diff(sampleId_minSnap));
disp("min-snap sparse GP jerk cost: "); disp(minSnapTraj_jerkCost);
