function [mu_query, sigma_query] = GPR_dense(x_query, x_data, y_data, eps_sigma, kernel_fcn)
%GPR_DENSE GP regression by dense kernel/covariance matrix 
%   GP form: Y_data = f(x_data) + eps, eps~N(0, eps_sigma*I)
%   Notes*: suppose Y reduced mean value, E[f(x)] = 0
%   inputs:
%       x_query: K X dim, query points 
%       x_data: N X dim, state data points 
%       y_data: N X dim, observation data points index
%       eps_sigma: noise term
%       kernel_fcn: kernel function handle -- kernel_fcn(x_data1, x_data2)

N = size(x_data, 1);
K_data = kernel_fcn(x_data, x_data) + eps_sigma*eye(N);
K_star = kernel_fcn(x_query, x_data);
K_query = kernel_fcn(x_query, x_query);

L = chol(K_data, "lower");
v = L\K_star;

mu_query = v' * (L\y_data);
sigma_query = K_query - v'*v;

end

