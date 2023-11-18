%test1: test state interpolation accuracy w.r.t time
clc; close all; clear;

%% Min-Acc 1-D tests
Qc = 0.01;
kappa0 = 0.01*eye(2);

statesSupp = [0, -10; 60, -300]; %random
t_supp = [0, 0.1];

testAccTraj = GPTrajSparseMinAcc(statesSupp, t_supp, Qc, statesSupp(1, :), kappa0);
testAccTraj.updateLifted();

figPos = [500, 100, 800, 600];

%test support Id range 1
figure("Position", figPos); 
subplot(2,1,1);
[sampStates, sampId] = plotSparseGP_1D(gca, testAccTraj, 1, 1000);
title("pos");

subplot(2,1,2);
plotSparseGP_1D(gca, testAccTraj, 2, 1000);
title("vel");

testPos = cumtrapz(sampId, sampStates(:, 2)) + sampStates(1, 1);
disp("min-Acc range 1 interp pos accuracy:");
error = max(abs(testPos - sampStates(:, 1)));
disp(error);

%test support Id range 2
testAccTraj.supportId = [0, 1];
testAccTraj.updateLifted();

figure("Position", figPos); 
subplot(2,1,1);
[sampStates, sampId] = plotSparseGP_1D(gca, testAccTraj, 1, 1000);
title("pos");

subplot(2,1,2);
plotSparseGP_1D(gca, testAccTraj, 2, 1000);
title("vel");

testPos = cumtrapz(sampId, sampStates(:, 2)) + sampStates(1, 1);
disp("min-Acc range 2 interp pos accuracy:");
error = max(abs(testPos - sampStates(:, 1)));
disp(error);

%test support Id range 3
testAccTraj.supportId = [0, 10];
testAccTraj.updateLifted();

figure("Position", figPos); 
subplot(2,1,1);
[sampStates, sampId] = plotSparseGP_1D(gca, testAccTraj, 1, 1000);
title("pos");

subplot(2,1,2);
plotSparseGP_1D(gca, testAccTraj, 2, 1000);
title("vel");

testPos = cumtrapz(sampId, sampStates(:, 2)) + sampStates(1, 1);
disp("min-Acc range 3 interp pos accuracy:");
error = max(abs(testPos - sampStates(:, 1)));
disp(error);

%% Min-Jerk 1-D tests
Qc = 0.01;
kappa0 = 0.01*eye(3);

statesSupp = [0, 10, 0; 30, 30, 0]; %random
t_supp = [0, 0.1];

testJerkTraj = GPTrajSparseMinJerk(statesSupp, t_supp, Qc, statesSupp(1, :), kappa0);
testJerkTraj.updateLifted();

figPos = [500, 100, 800, 900];

%test support Id range 1
figure("Position", figPos); 
subplot(3,1,1);
[sampStates, sampId] = plotSparseGP_1D(gca, testJerkTraj, 1, 1000);
title("pos");

subplot(3,1,2);
plotSparseGP_1D(gca, testJerkTraj, 2, 1000);
title("vel");

subplot(3,1,3);
plotSparseGP_1D(gca, testJerkTraj, 3, 1000);
title("acc");

testPos = cumtrapz(sampId, sampStates(:, 2)) + sampStates(1, 1);
disp("min-Jerk range 1 interp pos accuracy:");
error = max(abs(testPos - sampStates(:, 1)));
disp(error);

%test support Id range 2
testJerkTraj.supportId = [0, 1];
testJerkTraj.updateLifted();

figure("Position", figPos); 
subplot(3,1,1);
[sampStates, sampId] = plotSparseGP_1D(gca, testJerkTraj, 1, 1000);
title("pos");

subplot(3,1,2);
plotSparseGP_1D(gca, testJerkTraj, 2, 1000);
title("vel");

subplot(3,1,3);
plotSparseGP_1D(gca, testJerkTraj, 3, 1000);
title("acc");

testPos = cumtrapz(sampId, sampStates(:, 2)) + sampStates(1, 1);
disp("min-Jerk range 2 interp pos accuracy:");
error = max(abs(testPos - sampStates(:, 1)));
disp(error);

%test support Id range 3
testJerkTraj.supportId = [0, 10];
testJerkTraj.updateLifted();

figure("Position", figPos); 
subplot(3,1,1);
[sampStates, sampId] = plotSparseGP_1D(gca, testJerkTraj, 1, 1000);
title("pos");

subplot(3,1,2);
plotSparseGP_1D(gca, testJerkTraj, 2, 1000);
title("vel");

subplot(3,1,3);
plotSparseGP_1D(gca, testJerkTraj, 3, 1000);
title("acc");

testPos = cumtrapz(sampId, sampStates(:, 2)) + sampStates(1, 1);
disp("min-Jerk range 3 interp pos accuracy:");
error = max(abs(testPos - sampStates(:, 1)));
disp(error);

%% Min-Snap 1-D tests
Qc = 0.01;
kappa0 = 0.01*eye(4);

statesSupp = [0, 10, 100, 10; 30, 30, 10, -10]; %random
t_supp = [0, 0.1];

testSnapTraj = GPTrajSparseMinSnap(statesSupp, t_supp, Qc, statesSupp(1, :), kappa0);
testSnapTraj.updateLifted();

figPos = [500, 100, 1000, 800];

%test support Id range 1
figure("Position", figPos); 
subplot(2,2,1);
[sampStates, sampId] = plotSparseGP_1D(gca, testSnapTraj, 1, 1000);
title("pos");

subplot(2,2,2);
plotSparseGP_1D(gca, testSnapTraj, 2, 1000);
title("vel");

subplot(2,2,3);
plotSparseGP_1D(gca, testSnapTraj, 3, 1000);
title("acc");

subplot(2,2,4);
plotSparseGP_1D(gca, testSnapTraj, 4, 1000);
title("jerk");

testPos = cumtrapz(sampId, sampStates(:, 2)) + sampStates(1, 1);
disp("min-Snap range 1 interp pos accuracy:");
error = max(abs(testPos - sampStates(:, 1)));
disp(error);

%test support Id range 2
testSnapTraj.supportId = [0, 1];
testSnapTraj.updateLifted();

figure("Position", figPos); 
subplot(2,2,1);
[sampStates, sampId] = plotSparseGP_1D(gca, testSnapTraj, 1, 1000);
title("pos");

subplot(2,2,2);
plotSparseGP_1D(gca, testSnapTraj, 2, 1000);
title("vel");

subplot(2,2,3);
plotSparseGP_1D(gca, testSnapTraj, 3, 1000);
title("acc");

subplot(2,2,4);
plotSparseGP_1D(gca, testSnapTraj, 4, 1000);
title("jerk");

testPos = cumtrapz(sampId, sampStates(:, 2)) + sampStates(1, 1);
disp("min-Snap range 2 interp pos accuracy:");
error = max(abs(testPos - sampStates(:, 1)));
disp(error);

%test support Id range 3
testSnapTraj.supportId = [0, 10];
testSnapTraj.updateLifted();

figure("Position", figPos); 
subplot(2,2,1);
[sampStates, sampId] = plotSparseGP_1D(gca, testSnapTraj, 1, 1000);
title("pos");

subplot(2,2,2);
plotSparseGP_1D(gca, testSnapTraj, 2, 1000);
title("vel");

subplot(2,2,3);
plotSparseGP_1D(gca, testSnapTraj, 3, 1000);
title("acc");

subplot(2,2,4);
plotSparseGP_1D(gca, testSnapTraj, 4, 1000);
title("jerk");

testPos = cumtrapz(sampId, sampStates(:, 2)) + sampStates(1, 1);
disp("min-Snap range 3 interp pos accuracy:");
error = max(abs(testPos - sampStates(:, 1)));
disp(error);