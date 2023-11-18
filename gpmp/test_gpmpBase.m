%test basic functions in gpmpBase (for debugging)
clc; clear; close;

%% Setup test
% construct anchor states using cubic-spline
t_wp = [0, 5, 7, 8, 10, 15, 18];
q_wp = [3; -2; -5; 0; 6; 12; 8];
init_vel = 0;
end_vel = 0;
numSample = 10;

[q_cubic, dq_cubic, ddq_cubic, dddq_cubic, t_sample, ~]=cubicSplineTraj_wp(q_wp, t_wp, init_vel, end_vel, numSample);

% init prior params
Q_c = 1e-2; %noise variance at each dimension
kappa0 = diag([1e-1, 1e-1, 1e-1]); %init state variance
supportPts_minjerk = [q_cubic, dq_cubic, ddq_cubic];

minJerkTraj = GPTrajSparseMinJerk(supportPts_minjerk, t_sample, Q_c, supportPts_minjerk(1, :)', kappa0);

% init gpmp child class
options = gpmp_optimset_init();
varDim = 1; %1D
order = 3; %pos~acc
testObj = gpmpUncon2D(minJerkTraj, varDim, order, 0.1, 0.1, 0.1, options);

%% Test up sampling
%upsample by Matrix
testObj.upSampleByNums(ones(1,numSample-1) * 5);

supportPtsVec = testObj.getStateVecArrange(minJerkTraj.supportPts);

sampleStatesVec = testObj.getSampleStates(supportPtsVec);
sampleStates = testObj.getStateRowArrange(sampleStatesVec);

sampleTau = testObj.tau_samp;

%upsample by query
querySampleStates = zeros(size(sampleStates));
for i=1 : length(sampleTau)
    querySampleStates(i, :) = testObj.gpTrajSparse.query(sampleTau(i));
end

%compare
max(abs(sampleStates- querySampleStates), [], "all")
