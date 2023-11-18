%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function value = value_arclength_fn_wp( xi, c_fn_wp, p_start )
%VALUE_ARCLENGTH_FN_WP Calculate value of arc length cost function when
%given cost function handle
%   xi: discretized waypoint trajectory [Nxd]
%   c_fn: a function that maps a [1xd] point to a scalar
%   p_start: (optional) start point

if (nargin <= 2)
    p_start = xi(1,:);
end

c_fn_wpset = @(xi) stacked_fn( xi, c_fn_wp );
value = value_arclength_fn( xi, c_fn_wpset, p_start );

end

