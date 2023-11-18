%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function grad = gradient_arclength_fn_wp( xi, c_fn_wp, grad_c_fn_wp, p_start )
%GRADIENT_ARCLENGTH_FN_WP Calculate grad of arc length cost function when
%given cost function and gradient function handle
%   xi: discretized waypoint trajectory [Nxd]
%   c_fn_wp: a function that maps a [1xd] point to a scalar
%   grad_c_fn_wp: a function that maps a [1xd] point to a [1xd] gradient
%   p_start: (optional) start point

if (nargin <= 3)
    p_start = xi(1,:);
end

c_fn_wpset = @(xi) stacked_fn( xi, c_fn_wp );
grad_c_fn_wpset = @(xi) stacked_fn( xi, grad_c_fn_wp );
grad = gradient_arclength_fn( xi, c_fn_wpset, grad_c_fn_wpset, p_start);

end

