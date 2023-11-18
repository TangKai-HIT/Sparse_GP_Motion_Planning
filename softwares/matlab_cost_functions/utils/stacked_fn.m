%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function val = stacked_fn( xi, fn )
%STACKED_FN Call a function for every point in xi and stack the output
%   xi: Nxd set of waypoints, N is number of waypoints
%   fn: Function to be evaluated on each waypoint
%   val: stacked output of each function call

val = [];
for t = 1:size(xi,1)
    val = [val; fn(xi(t,:))];
end

end

