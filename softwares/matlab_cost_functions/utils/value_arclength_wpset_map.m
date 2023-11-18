%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function value = value_arclength_wpset_map( wpset, map, p_start )
%VALUE_ARCLENGTH_WPSET_MAP Arclength integration to get value over a wpset
%   wpset: Nx2 waypoint set, N is number of waypoints
%   map: default map struct
%   p_start: (optional) start coordinate to anchor trajectory
%   value: scalar arclength value along trajectory


if (nargin <= 2)
    p_start = wpset(1,:);
end

value_set = value_wpset_map( wpset, map );
value = sum(sqrt(sum((diff([p_start; wpset])).^2,2)) .* value_set);

end
