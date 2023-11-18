%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function cost = costfn_geometric_value_wpset_map( wpset, map, infinity )
%COSTFN_GEOMETRIC_VALUE_WPSET_MAP Get value of a wp traj for geometric cost
%function, i.e, length if the line is collision free, else 
%   wpset: Nx2 waypoint set, N is number of waypoints
%   map: default map struct
%   infinity: value if in collision (defaults to inf)
%   cost: cost of the line

if (nargin == 3)
    infinity = inf;
end

if (~check_coll_wpset_map( wpset, map ))
    cost = infinity;
else
    cost = costfn_length_value_wpset( wpset );
end

end

