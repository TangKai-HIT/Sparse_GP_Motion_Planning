function K = RBF_kernel(x_data1, x_data2, sigma_f, len)
%RBF_KERNEL RBF kernel
%   inputs:
%       x_data1, x_data2: k1 X dim, k2 X dim
%       sigma_f, len: hyperparams

M = pdist2(x_data1, x_data2, "squaredeuclidean");

K = sigma_f * exp(- M ./ (2*len^2));

end

