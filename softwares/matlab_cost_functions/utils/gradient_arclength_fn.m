%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function grad = gradient_arclength_fn( xi, c_fn, grad_c_fn, p_start )
%GRADIENT_ARCLENGTH_FN Calculate gradient of arc length cost function when
%given cost function and gradient function handles
%   xi: discretized waypoint trajectory [Nxd]
%   c_fn: a function that maps a [Nxd] wpset to a [Nx1] value
%   grad_c_fn: a function that maps a [Nxd] wpset to a [Nxd] gradient value
%   p_start: (optional) start point

if (nargin <= 3)
    p_start = xi(1,:);
end

n = size(xi,1) + 1;
d = size(xi,2);
xi_d = n*diff([p_start; xi]);
xi_d_norm = normr(xi_d);
xi_dd = n*([zeros(1,d); diff(xi_d)]);
kappa = repmat((1./sum((xi_d).^2,2)), [1 d]).*(xi_dd - xi_d_norm.*repmat(sum(xi_dd.*xi_d_norm, 2), [1 d]));

c = c_fn(xi);
delta_c = grad_c_fn(xi);
grad = (1/n).*repmat(sqrt(sum((xi_d).^2,2)), [1 d]) .* ( delta_c - xi_d_norm.*repmat(sum(delta_c.*xi_d_norm,2), [1 d])  - repmat(c, [1 d]).*kappa);
grad(isnan(grad)) = 0;

end

