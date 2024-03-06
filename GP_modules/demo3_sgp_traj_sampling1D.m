%demo3: sample sparse 1-D GP trajectory between given boundary states
close all; clear; clc;

%% Define support states
statesSupp = [0, 10, 0; 30, 30, 0]; %random
t_supp = [0, 10];

%% Init prior params
dim = 3;
Q_c = 1e2; %noise variance at each dimension
kappa0 = diag([1e-1, 1e-1, 1e-1]); %init state variance

minJerkTraj = GPTrajSparseMinJerk(statesSupp, t_supp, Q_c, statesSupp(1, :), kappa0);

minJerkTraj.updateLifted();

%% Get distribution on sample time index
%samples 1
numStates1 = 3;
sampTimeIdAll1 = linspace(t_supp(1), t_supp(end), numStates1+2);
sampTimeId1 = sampTimeIdAll1(2:end-1);

[~, ~, ~, sampDist1_mu, sampDist1_kappa, sampDist1_kappa_inv] = minJerkTraj.upSampleInterval(1, 2, sampTimeId1);
sampDist1_mu_row = getStateRowArrange(sampDist1_mu, 3);
sampDist1_kappa = (sampDist1_kappa + sampDist1_kappa')/2;

minJerkSampTraj1 = GPTrajSparseMinJerk([statesSupp(1,:); sampDist1_mu_row; statesSupp(end,:)], sampTimeIdAll1, Q_c, statesSupp(1, :), kappa0);
minJerkSampTraj1.updateLifted();

%samples 2
numStates2 = 5;
sampTimeIdAll2 = linspace(t_supp(1), t_supp(end), numStates2+2);
sampTimeId2 = sampTimeIdAll2(2:end-1);

[~, ~, ~, sampDist2_mu, sampDist2_kappa, sampDist2_kappa_inv] = minJerkTraj.upSampleInterval(1, 2, sampTimeId2);
sampDist2_mu_row = getStateRowArrange(sampDist2_mu, 3);
sampDist2_kappa = (sampDist2_kappa + sampDist2_kappa')/2;

minJerkSampTraj2 = GPTrajSparseMinJerk([statesSupp(1,:); sampDist2_mu_row; statesSupp(end,:)], sampTimeIdAll2, Q_c, statesSupp(1, :), kappa0);
minJerkSampTraj2.updateLifted();

%% Sample trajs
rng(3);

%sample 1
numSamples1 = 10;
sampleStates1_row = mvnrnd(sampDist1_mu', sampDist1_kappa, numSamples1);

%sample 2
numSamples2 = 10;
sampleStates2_row = mvnrnd(sampDist2_mu', sampDist2_kappa, numSamples2);

%% Plot sampled SGP trajectory 1
N = 500;
%get states on sampled trajectory 1
sampleStates1 = zeros(N, dim, numSamples1);
for i=1:numSamples1
    curSampStates1 = getStateRowArrange(sampleStates1_row(i, :)', 3);
    minJerkSampTraj1.supportPts = [statesSupp(1,:); curSampStates1; statesSupp(end,:)];
    minJerkSampTraj1.updateLifted();
    [sampleStates1(:, :, i), sampleId1] = sampleStatesHelper(minJerkSampTraj1, N);
end

%set color map
cmap1 = cool(numSamples1);

%set var scale 
varScale = 2;

%plots
figure('Position',[500,80, 1000,900]);
subplot(3,1,1);
plotSparseGP_1D(gca, minJerkTraj, 1, N, varScale);
title("position");
hold on;
for i=1:numSamples1
    plot(sampleId1, sampleStates1(:, 1, i), '-', LineWidth=0.8, Color=cmap1(i, :));
    scatter(sampTimeId1, sampleStates1_row(i, 1:dim:end), MarkerFaceColor=cmap1(i, :), LineWidth=0.4);
end
legend("conditioned variance","mean prediction", "support states", "Location","northeastoutside");

subplot(3,1,2);
plotSparseGP_1D(gca, minJerkTraj, 2, N, varScale);
title("velocity");
hold on;
for i=1:numSamples1
    plot(sampleId1, sampleStates1(:, 2, i), '-', LineWidth=0.8, Color=cmap1(i, :));
    scatter(sampTimeId1, sampleStates1_row(i, 2:dim:end), MarkerFaceColor=cmap1(i, :), LineWidth=0.4);
end
legend("conditioned variance","mean prediction", "support states", "Location","northeastoutside");

subplot(3,1,3);
plotSparseGP_1D(gca, minJerkTraj, 3, N, varScale);
title("acceleration");
hold on;
for i=1:numSamples1
    plot(sampleId1, sampleStates1(:, 3, i), '-', LineWidth=0.8, Color=cmap1(i, :));
    scatter(sampTimeId1, sampleStates1_row(i, 3:dim:end), MarkerFaceColor=cmap1(i, :), LineWidth=0.4);
end
legend("conditioned variance","mean prediction", "support states", "Location","northeastoutside");

sgtitle("SGP trajectory samples 1");

%% Plot sampled SGP trajectory 2
N = 500;
%get states on sampled trajectory 2
sampleStates2 = zeros(N, dim, numSamples2);
for i=1:numSamples2
    curSampStates2 = getStateRowArrange(sampleStates2_row(i, :)', 3);
    minJerkSampTraj2.supportPts = [statesSupp(1,:); curSampStates2; statesSupp(end,:)];
    minJerkSampTraj2.updateLifted();
    [sampleStates2(:, :, i), sampleId2] = sampleStatesHelper(minJerkSampTraj2, N);
end

%set color map
cmap2 = cool(numSamples2);

%set var scale 
varScale = 2;

%plots
figure('Position',[500,80, 1000,900]);
subplot(3,1,1);
plotSparseGP_1D(gca, minJerkTraj, 1, N, varScale);
title("position");
hold on;
for i=1:numSamples2
    plot(sampleId2, sampleStates2(:, 1, i), '-', LineWidth=0.8, Color=cmap2(i, :));
    scatter(sampTimeId2, sampleStates2_row(i, 1:dim:end), MarkerFaceColor=cmap2(i, :), LineWidth=0.4);
end
legend("conditioned variance","mean prediction", "support states", "Location","northeastoutside");

subplot(3,1,2);
plotSparseGP_1D(gca, minJerkTraj, 2, N, varScale);
title("velocity");
hold on;
for i=1:numSamples2
    plot(sampleId2, sampleStates2(:, 2, i), '-', LineWidth=0.8, Color=cmap2(i, :));
    scatter(sampTimeId2, sampleStates2_row(i, 2:dim:end), MarkerFaceColor=cmap2(i, :), LineWidth=0.4);
end
legend("conditioned variance","mean prediction", "support states", "Location","northeastoutside");

subplot(3,1,3);
plotSparseGP_1D(gca, minJerkTraj, 3, N, varScale);
title("acceleration");
hold on;
for i=1:numSamples2
    plot(sampleId2, sampleStates2(:, 3, i), '-', LineWidth=0.8, Color=cmap2(i, :));
    scatter(sampTimeId2, sampleStates2_row(i, 3:dim:end), MarkerFaceColor=cmap2(i, :), LineWidth=0.4);
end
legend("conditioned variance","mean prediction", "support states", "Location","northeastoutside");

sgtitle("SGP trajectory samples 2");
