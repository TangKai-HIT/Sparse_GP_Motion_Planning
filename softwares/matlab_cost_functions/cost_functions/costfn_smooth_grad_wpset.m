%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function grad = costfn_smooth_grad_wpset( xi, A, b )
%GRADIENT_SMOOTH Calculate smoothness gradient
%   xi: discretized waypoint trajectory
%   A: smoothness quadratic matrix
%   b: smoothness linear matrix
%   cost: smoothness cost

% grad = size(xi,1)*(A*xi+b); %original
grad = (size(xi,1)+1)*(A*xi+b); %modified
end

