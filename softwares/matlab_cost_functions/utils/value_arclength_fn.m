%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function value = value_arclength_fn( xi, c_fn, p_start )
%VALUE_ARCLENGTH_FN Calculate value of arc length cost function when
%given cost function handle
%   xi: discretized waypoint trajectory [Nxd]
%   c_fn: a function that maps a [Nxd] point to a [Nx1] vec
%   p_start: (optional) start point

if (nargin <= 2)
    p_start = xi(1,:);
end

value = sum(sqrt(sum((diff([p_start; xi])).^2,2)).* c_fn(xi));
end

