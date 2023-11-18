function [sampleStates, sampleId] = plotSparseGP_2D(ax, GPTrajSparse, dimx, dimy, N)
%PLOTSPARSEGP_2D 2-D trajectory plot util for 2-D GPTrajSparse objects
%   GPTrajSparse: derived-class from GPTrajSparseBase
%   dimx, dimy: x-axis dim index, y-axis dim index

sampleId = linspace(GPTrajSparse.supportId(1), GPTrajSparse.supportId(end), N);
sampleStates = zeros(N, GPTrajSparse.dim);
% sampleVars = zeros(GPTrajSparse.dim, GPTrajSparse.dim, N);

for i = 1:N
    sampleStates(i, :) = GPTrajSparse.query(sampleId(i));
end

hold(ax, "on");

% plot expectation trajectory
plot(ax, sampleStates(:, dimx), sampleStates(:, dimy), '-r', 'LineWidth', 1.5);
scatter(ax, GPTrajSparse.supportPts(:, dimx), GPTrajSparse.supportPts(:, dimy), 'black', 'filled' ,'o', SizeData=30);

end

