%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function [ cost_map ] = create_costmap_sqdist( map, epsilon )
%CREATE_COSTMAP_SQDIST Create a cost map that computes squared distance to
%obstacles (linear distance when inside)
%   map: environment map struct
%   epsilon: horizon upto which distance is computed
%   cost_map: cost_map where each coordinate maps to the cost

Dint = double(-bwdist(map.table));
Dext = double(bwdist(1-map.table));

Cint = -map.resolution*Dint;
Cext = (1.0/(2.0*epsilon))*((min(Dext*map.resolution-epsilon, 0)).^2);

cost_map = map;
cost_map.table = Cint + Cext;

end

