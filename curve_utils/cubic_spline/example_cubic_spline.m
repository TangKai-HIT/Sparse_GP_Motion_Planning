% A simple cubic spline example
% addpath(genpath(pwd));
%% Init 
t_wp = [0, 5, 7, 8, 10, 15, 18];
q_wp = [3; -2; -5; 0; 6; 12; 8];
init_vel = 2;
end_vel = -3;
numSample = 100;

%% Compute trajectory
[q, dq, ddq, dddq, t_sample, splineCoeff]=cubicSplineTraj_wp(q_wp, t_wp, init_vel, end_vel, numSample);

%% Plot
plotViaPtSpline1D(t_sample, q, dq, ddq, dddq, t_wp, q_wp);
