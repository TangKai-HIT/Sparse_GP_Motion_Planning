%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function [ cost ] = costfn_geometric_value_line_map( p_start, p_end, map, infinity )
%COSTFN_GEOMETRIC_VALUE_LINE_MAP Get value of a line for geometric cost
%function, i.e, length if the line is collision free, else infinity
%   p_start: 1xd point
%   p_end: 1xd point
%   map: map struct
%   infinity: value if in collision (defaults to inf)
%   cost: cost of the line

if (nargin == 3)
    infinity = inf;
end

if (~check_coll_line_map(p_start, p_end, map))
    cost = infinity;
else
    cost = norm(p_end - p_start);
end

end

