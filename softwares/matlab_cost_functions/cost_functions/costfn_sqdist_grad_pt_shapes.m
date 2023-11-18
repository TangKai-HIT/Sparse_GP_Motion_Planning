%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function grad = costfn_sqdist_grad_pt_shapes( pt, shapes_array, epsilon )
%COSTFN_SQDIST_GRAD_PT_SHAPES Get grad of a point
%   pt: 1xd point
%   shapes_array: array of shapes
%   epsilon: horizon of expansion
%   c: cost at the point
%   c_grad: gradient of cost

[~,grad] = costfn_sqdist_pt_shapes( pt, shapes_array, epsilon );

end
