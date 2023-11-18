%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function grad = gradient_arclength_wpset_map( wpset, cost_map, cost_map_x, cost_map_y, p_start )
%GRADIENT_ARCLENGTH_WPSET_MAP Get gradient of arc length cost functions
%   wpset: Nx2 waypoint set, N is number of waypoints
%   cost_map: default map struct
%   cost_map_x: x gradient of cost_map
%   cost_map_y: y gradient of cost_map
%   p_start: (optional) start coordinate to anchor trajectory
%   grad: gradient of trajectory

if (nargin <= 4)
    p_start = wpset(1,:);
end

xi = wpset;

c = value_wpset_map( xi, cost_map );
delta_c = [value_wpset_map( xi, cost_map_x ) value_wpset_map( xi, cost_map_y )];

n = size(xi,1) + 1;
xi_d = n*diff([p_start; xi]);
xi_d_norm = normr(xi_d);
xi_dd = n*([0 0; diff(xi_d)]);

kappa = repmat((1./sum((xi_d).^2,2)), [1 2]).*(xi_dd - xi_d_norm.*repmat(sum(xi_dd.*xi_d_norm, 2), [1 2]));

grad = (1/n).*repmat(sqrt(sum((xi_d).^2,2)), [1 2]) .* ( delta_c - xi_d_norm.*repmat(sum(delta_c.*xi_d_norm,2), [1 2])  - repmat(c, [1 2]).*kappa);
grad(isnan(grad)) = 0;
end

