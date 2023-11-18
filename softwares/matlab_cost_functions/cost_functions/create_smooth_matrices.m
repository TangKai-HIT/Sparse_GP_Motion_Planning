%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function [ A, b, c, K, e ] = create_smooth_matrices( p_start, p_goal, n )
%CREATE_SMOOTH_MATRICES Create matrices to compute smooth cost function
%   p_start: anchor start point
%   p_goal: anchor end point
%   n: number of points of trajectory
%   A: smoothness quadratic matrix
%   b: smoothness linear matrix
%   c: smoothness constant matrix
%   K: differentiation matrix
%   e: error matrix

%% Set up smoothness matrices
% f = 0.5*||Kxi + e||^2
K = zeros(n+1, n);
fd = [-1 1 zeros(1, n-2)];
for i = 1:(n-1)
    K(i+1,:) = circshift(fd, [0 (i-1)]);
end
K(1,1) = 1;
K(n+1,n) = -1;

e = [-p_start; zeros(n-1, length(p_start)); p_goal];

A = K'*K;
b = K'*e;
c = 0.5*(e'*e);

end

