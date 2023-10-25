function [sampleStates, sampleId] = plotSparseGP_1D(ax, GPTrajSparse, dim, N, varScale)
%PLOTSPARSEGP_1D 1-D trajectory plot util for GPTrajSparse objects
%   GPTrajSparse: derived-class from GPTrajSparseBase
%   varScale: scale scalar of state variance plot

if ~exist("varScale", "var")
    varScale = 2E1;
end

sampleId = linspace(GPTrajSparse.supportId(1), GPTrajSparse.supportId(end), N);
sampleStates = zeros(N, GPTrajSparse.dim);
sampleVars = zeros(GPTrajSparse.dim, GPTrajSparse.dim, N);

for i = 1:N
    [sampleStates(i, :), sampleVars(:, :, i)] = GPTrajSparse.query(sampleId(i));
end

hold(ax, "on");
% plot variance patch
patch(ax, [sampleId, sampleId(end:-1:1)], ...
		[sampleStates(:, dim)+squeeze(abs(sampleVars(dim, dim, :)).^.5)*varScale; sampleStates(end:-1:1, dim)-squeeze(abs(sampleVars(dim, dim, end:-1:1)).^.5)*varScale], ...
		[1 .8 .8],'edgecolor','none'); %added abs to avoid negative value around 0 (numerical problem) 

% plot expectation trajectory
plot(ax, sampleId, sampleStates(:, dim), '-r', 'LineWidth', 1.5);
scatter(ax, GPTrajSparse.supportId, GPTrajSparse.supportPts(:, dim), 'black', 'filled' ,'o', SizeData=30);

min_y = min(sampleStates(:, dim));
max_y = max(sampleStates(:, dim));
mid_y = (min_y + max_y)/2;
gap = max_y - min_y;
ylim([mid_y-3*gap/2, mid_y+3*gap/2]);
grid on;

end

