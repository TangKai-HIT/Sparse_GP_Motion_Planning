function [sampleStates, sampleId] = sampleStatesHelper(GPTrajSparse, N)
%SAMPLESTATESHELPER sample states in sparse gp traj
%   此处显示详细说明

sampleId = linspace(GPTrajSparse.supportId(1), GPTrajSparse.supportId(end), N);
sampleStates = zeros(N, GPTrajSparse.dim);

for i = 1:N
    sampleStates(i, :) = GPTrajSparse.query(sampleId(i));
end

end

