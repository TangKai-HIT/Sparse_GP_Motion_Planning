%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function [ cost_map_x, cost_map_y, cost_map_xx, cost_map_yy, cost_map_xy ] = get_cost_map_derivatives( cost_map )
%GET_COST_MAP_DERIVATIVES Gets derivatives and double derivatives of cost
%map
%   cost_map: cost map to find derivative of
%   cost_map_x: x derivative
%   cost_map_y: y derivative
%   cost_map_xx: xx derivative
%   cost_map_yy: yy derivative
%   cost_map_xy: xy derivative

cost_map_x = cost_map;
cost_map_y = cost_map;
cost_map_xx = cost_map;
cost_map_yy = cost_map;
cost_map_xy = cost_map;

[Cy,Cx] = gradient(cost_map.table/cost_map.resolution);
[Cyx, Cxx] = gradient(Cx/cost_map.resolution);
[Cyy, Cxy] = gradient(Cy/cost_map.resolution);

cost_map_x.table = Cx;
cost_map_y.table = Cy;
cost_map_xx.table = Cxx;
cost_map_yy.table = Cyy;
cost_map_xy.table = 0.5*(Cxy + Cyx); %forcing symmetry
end

