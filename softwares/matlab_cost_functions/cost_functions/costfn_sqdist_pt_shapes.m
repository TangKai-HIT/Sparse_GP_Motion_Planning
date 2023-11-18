%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function [ c, c_grad ] = costfn_sqdist_pt_shapes( pt, shapes_array, epsilon )
%COSTFN_SQDIST_PT_SHAPES Get value and gradient of a point
%   pt: 1xd point
%   shapes_array: array of shapes
%   epsilon: horizon of expansion
%   c: cost at the point
%   c_grad: gradient of cost

[~, d, grad] = shapes_point_check( pt, shapes_array );
if (d < 0)
    c = -d + 0.5*epsilon;
    c_grad = -grad;
elseif (d < epsilon)
    c = (1/(2*epsilon))*(d - epsilon)^2;
    c_grad = (1/epsilon)*(d - epsilon)*grad;
else
    c = 0;
    c_grad = zeros(size(pt));
end

end

