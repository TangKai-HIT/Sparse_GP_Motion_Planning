%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function cost = costfn_smooth_value_wpset( xi, A, b, c )
%COSTFN_SMOOTH_VALUE_WPSET Calculate smoothness of a trajectory
%   xi: discretized waypoint trajectory
%   A: smoothness quadratic matrix
%   b: smoothness linear matrix
%   c: smoothness constant matrix
%   cost: smoothness cost

cost = (size(xi,1)+1)*trace(0.5*xi'*A*xi + xi'*b + c);

end

