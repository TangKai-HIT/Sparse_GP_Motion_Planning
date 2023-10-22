function sparse_mat = getTridiag(tri_mat)
%GETTRIDIAG extract sparse elements in a tridiagonal matrix to a N x 3 tuple array
%   Inputs:
%       tri_mat: N x N, tridiagonal matrix, N>=2
%   outputs:
%       sparse_mat: N x 3 matrix with tridiagonal elements stored in each row -- [a_n, b_n, c_n]

N = size(tri_mat, 1);
sparse_mat = zeros(N, 3);

sparse_mat(1, 2:3) = tri_mat(1, 1:2);
sparse_mat(N, 1:2) = tri_mat(N, N-1:N);

for i=2:N-1
    sparse_mat(i, :) = tri_mat(i, i-1:i+1);
end

