%test2: compare error and acc/jerk costs with polynomial trajectories
clc; close all; clear;

%% Min-Acc 1-D tests (compare with cubic polynomial)
Qc = 0.01;
kappa0 = 0.01*eye(2);

statesSupp = [0, -10; 60, -300]; %random
t_supp = [0, 1];

testAccTraj = GPTrajSparseMinAcc(statesSupp, t_supp, Qc, statesSupp(1, :), kappa0);
testAccTraj.updateLifted();

figPos = [500, 50, 800, 600];

%test support Id range 1
figure("Position", figPos); 
subplot(2,1,1);
[sampStates, t_samp_quin] = plotSparseGP_1D(gca, testAccTraj, 1, 1000);
[q,qd,~] = cubicpolytraj(statesSupp(:,1)', testAccTraj.supportId, t_samp_quin, ...
                            "VelocityBoundaryCondition", statesSupp(:, 2)');
hold on;
plot(t_samp_quin, q, '--b');
title("pos");

subplot(2,1,2);
plotSparseGP_1D(gca, testAccTraj, 2, 1000);
hold on;
plot(t_samp_quin, qd, '--b');
title("vel");

disp("min-Acc vs cubic-poly range 1 pos error:");
error = max(abs(q' - sampStates(:, 1)));
disp(error);

%test support Id range 2
testAccTraj.supportId = [0, 10];
testAccTraj.updateLifted();

figure("Position", figPos); 
subplot(2,1,1);
[sampStates, t_samp_quin] = plotSparseGP_1D(gca, testAccTraj, 1, 1000);
[q,qd,~] = cubicpolytraj(statesSupp(:,1)', testAccTraj.supportId, t_samp_quin, ...
                            "VelocityBoundaryCondition", statesSupp(:, 2)');
hold on;
plot(t_samp_quin, q, '--b');
title("pos");

subplot(2,1,2);
plotSparseGP_1D(gca, testAccTraj, 2, 1000);
hold on;
plot(t_samp_quin, qd, '--b');
title("vel");

disp("min-Acc  vs cubic-poly range 2 interp pos accuracy:");
error = max(abs(q' - sampStates(:, 1)));
disp(error);

%% Min-Jerk 1-D tests (compare with quintic polynomial)
Qc = 0.01;
kappa0 = 0.01*eye(3);

statesSupp = [0, 10, 0; 30, 30, 0]; 
t_supp = [0, 1];

testJerkTraj = GPTrajSparseMinJerk(statesSupp, t_supp, Qc, statesSupp(1, :), kappa0);
testJerkTraj.updateLifted();

figPos = [500, 50, 800, 900];

%test support Id range 1
figure("Position", figPos); 
subplot(3,1,1);
[sampStates, sampId] = plotSparseGP_1D(gca, testJerkTraj, 1, 1000);
[q,qd,qdd] = quinticpolytraj(statesSupp(:,1)', testJerkTraj.supportId, sampId, ...
                            "VelocityBoundaryCondition", statesSupp(:, 2)', "AccelerationBoundaryCondition", statesSupp(:, 3)');
hold on;
plot(sampId, q, '--b');
title("pos");

subplot(3,1,2);
plotSparseGP_1D(gca, testJerkTraj, 2, 1000);
hold on;
plot(sampId, qd, '--b');
title("vel");

subplot(3,1,3);
plotSparseGP_1D(gca, testJerkTraj, 3, 1000);
hold on;
plot(sampId, qdd, '--b');
title("acc");

disp("min-Jerk vs quintic-poly range 1 interp pos accuracy:");
error = max(abs(q' - sampStates(:, 1)));
disp(error);

%test support Id range 2
testJerkTraj.supportId = [0, 10];
testJerkTraj.updateLifted();

figure("Position", figPos); 
subplot(3,1,1);
[sampStates, sampId] = plotSparseGP_1D(gca, testJerkTraj, 1, 1000);
[q,qd,qdd] = quinticpolytraj(statesSupp(:,1)', testJerkTraj.supportId, sampId, ...
                            "VelocityBoundaryCondition", statesSupp(:, 2)', "AccelerationBoundaryCondition", statesSupp(:, 3)');
hold on;
plot(sampId, q, '--b');
title("pos");

subplot(3,1,2);
plotSparseGP_1D(gca, testJerkTraj, 2, 1000);
hold on;
plot(sampId, qd, '--b');
title("vel");

subplot(3,1,3);
plotSparseGP_1D(gca, testJerkTraj, 3, 1000);
hold on;
plot(sampId, qdd, '--b');
title("acc");

disp("min-Jerk vs quintic-poly range 2 interp pos accuracy:");
error = max(abs(q' - sampStates(:, 1)));
disp(error);
