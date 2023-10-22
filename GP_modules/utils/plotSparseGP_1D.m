function [sampleStates, sampleId] = plotSparseGP_1D(ax, GPTrajSparse)
%PLOTSPARSEGP_1D 1-D trajectory plot util for GPTrajSparse objects
%   GPTrajSparse: derived-class from GPTrajSparseBase

N = 500;
sampleId = linspace(GPTrajSparse.supportId(1), GPTrajSparse.supportId(end), N);
sampleStates = zeros(N, GPTrajSparse.dim);
sampleVars = zeros(GPTrajSparse.dim, GPTrajSparse.dim, N);

for i = 1:N
    [sampleStates(i, :), sampleVars(:, :, i)] = GPTrajSparse.query(sampleId(i));
end

hold(ax, "on");
% plot variance patch
patch(ax, [sampleId, sampleId(end:-1:1)], ...
		[sampleStates(:, 1)+squeeze(abs(sampleVars(1, 1, :)).^.5)*2E1; sampleStates(end:-1:1, 1)-squeeze(abs(sampleVars(1, 1, end:-1:1)).^.5)*2E1], ...
		[1 .8 .8],'edgecolor','none'); %added abs to avoid negative value around 0 (numerical problem) 

% plot expectation trajectory
plot(ax, sampleId, sampleStates(:, 1), '-r', 'LineWidth', 1.5);
scatter(ax, GPTrajSparse.supportId, GPTrajSparse.supportPts(:, 1), 'black', 'filled' ,'o', SizeData=30);

min_y = min(sampleStates(:, 1));
max_y = max(sampleStates(:, 1));
mid_y = (min_y + max_y)/2;
gap = max_y - min_y;
ylim([min_y-gap/2, max_y+gap/2]);
grid on;

end

