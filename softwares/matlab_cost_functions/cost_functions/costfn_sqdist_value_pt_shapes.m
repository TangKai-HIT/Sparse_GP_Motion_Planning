%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function value = costfn_sqdist_value_pt_shapes( pt, shapes_array, epsilon )
%COSTFN_SQDIST_VALUE_PT_SHAPES Get value of a point
%   pt: 1xd point
%   shapes_array: array of shapes
%   epsilon: horizon of expansion
%   c: cost at the point
%   c_grad: gradient of cost

value = costfn_sqdist_pt_shapes( pt, shapes_array, epsilon );

end

