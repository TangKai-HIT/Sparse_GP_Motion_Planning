%demo2: construct sparse GP trajectory on SO(3)
close all; clear; clc;

%% Init test trajectory on SO(3) 
% construct anchor states
r0 = [1 0 0; 0 1 0; 0 0 1];
r1 = eul2rotm([0,pi/3,0], "ZXY"); % note: rotation on Z-axis is hard to show on the sphere
r2 = eul2rotm([0,pi/3,pi/2], "ZXY");

tSamples1 = linspace(0, 10, 250);
[R_sq1, omega_sq1, alpha_sq1] = rottraj(r0, r1, [0, 10], tSamples1);

tSamples2 = linspace(10, 20, 252);
tSamples2(1) = [];
[R_sq2, omega_sq2, alpha_sq2] = rottraj(r1, r2, [10, 20], tSamples2);

tSamples = [tSamples1, tSamples2];
R_sq = cat(3, R_sq1, R_sq2);
omega_sq = [omega_sq1, omega_sq2];

%% Show test path on sphere
figure('Position',[100,100, 800,800], "Name", "test path");
plotTestSamp_id = 1: 10: size(R_sq, 3);

% plot SO3 frames 
testSO3_path = so3(R_sq(:, :, plotTestSamp_id));
translations = zeros(length(plotTestSamp_id), 3);
for i=1:length(plotTestSamp_id)
    translations(i, :) = R_sq(:, :, plotTestSamp_id(i)) * [0;0;1];
end

testPath_ax1 = plotTransforms(translations,testSO3_path', 'FrameSize', 0.4);

hold on; axis off; grid off; 
% plot scatter path
scatter3(testPath_ax1, translations(:,1), translations(:,2), translations(:,3), "red", "filled", "o", "LineWidth", 0.5);

% plot sphere
[X,Y,Z] = sphere(30);
mesh(testPath_ax1, X,Y,Z, "EdgeColor", "black");
view(60,30); axis equal; axis vis3d;  rotate3d on;

title("test path");

%% Show test path on line
figure('Position',[100,100, 1100,800], "Name", "test path on line");
translations = zeros(length(plotTestSamp_id), 3);
translations(:, 2) = linspace(0, 30, length(plotTestSamp_id));
testPath_ax2 = plotTransforms(translations,testSO3_path', 'FrameSize', 1);

hold on; 
% plot scatter path
scatter3(testPath_ax2, translations(:,1), translations(:,2), translations(:,3), "red", "filled", "o", "LineWidth", 0.5);

title("test path on line");
axis equal; rotate3d on; view(20, 15);

%% Init min-acc sparse GP trajectory on SO(3)
sampleStep = 100;
t_support_id = 1:sampleStep:length(tSamples);
t_support = tSamples(t_support_id);

support_SO3 = R_sq(:, :, t_support_id);
support_so3 = omega_sq(:, t_support_id);

% init prior params
Q_c = 1e-2 * eye(3); %noise variance at each dimension
kappa0 = diag([1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1]); %init state variance

minAcc_SO3Traj = GPTrajSparseSO3MinAcc(support_SO3, support_so3, t_support, Q_c, kappa0);

minAcc_SO3Traj.updateLifted();

%% Plot min-Acc sparse GP trajectory SO(3)
sampId = linspace(minAcc_SO3Traj.supportId(1), minAcc_SO3Traj.supportId(end), 50);
sampSO3 = zeros(3,3,length(sampId));
sampOmega = zeros(3,length(sampId));

for i=1:length(sampId)
    [sampSO3(:, :, i), sampOmega(:, i)] = minAcc_SO3Traj.query(sampId(i));
end

figure('Position',[1000,100, 800,800], "Name", "min-Acc sparse GP mean");

% plot SO3 frames 
gpSO3_path = so3(sampSO3);
translations = zeros(length(sampId), 3);
for i=1:length(sampId)
    translations(i, :) = sampSO3(:, :, i) * [0;0;1];
end

sampPath_ax1 = plotTransforms(translations,gpSO3_path', 'FrameSize', 0.4);
hold on; axis off; grid off; 

% plot scatter path
scatter3(sampPath_ax1, translations(:,1), translations(:,2), translations(:,3), "red", "filled", "o", "LineWidth", 0.5);

% plot sphere
mesh(sampPath_ax1, X,Y,Z, "EdgeColor", "black");
view(60,30); axis equal; axis vis3d;  rotate3d on;

title("min-Acc sparse GP mean");

%% Show min-Acc sparse GP SO(3) on line
figure('Position',[100,100, 1100,800], "Name", "min-Acc sparse GP mean on line");
translations = zeros(length(sampId), 3);
translations(:, 2) = linspace(0, 30, length(sampId));
sampPath_ax2 = plotTransforms(translations,gpSO3_path', 'FrameSize', 1);

hold on; 
% plot scatter path
scatter3(sampPath_ax2, translations(:,1), translations(:,2), translations(:,3), "red", "filled", "o", "LineWidth", 0.5);

title("min-Acc sparse GP mean on line");
axis equal; rotate3d on; view(20, 15);