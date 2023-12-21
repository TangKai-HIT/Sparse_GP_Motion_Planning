function result = getDiffMatrix_quadInterp(points, order)
%GETDIFFMATRIX_QUADINTERP get differential matrix of unevenly sampled points using symmetric quadratic polynomial interpolation
%   Inputs:
%       points: N X 1 or 1 X N, sampled points
%       order: differential order (lower than 2nd order for quadratic interpolation)
%   Output:
%       result: N X N

N = length(points);
delta_s = diff(points, 1);
result = zeros(N, N);

for i = 1:N
    if i==1
        switch order
            case 1
                result(i, :) = [-1, 1, zeros(1, N-2)] / delta_s(1);
            case 2
                result(i, :) = [1/delta_s(1)^2,  -(1/(delta_s(1)*delta_s(2)) + 1/delta_s(1)^2), 1/(delta_s(1)*delta_s(2)), zeros(1, N-3)];
            otherwise
                disp('order <=2!');
                result = [];
                break;
        end

    elseif i==N
        switch order
            case 1
                result(i, :) = [zeros(1, N-2), -1, 1] / delta_s(end);
            case 2
                result(i, :) = [zeros(1, N-3), 1/(delta_s(end)*delta_s(end-1)),  -(1/(delta_s(end-1)*delta_s(end)) + 1/delta_s(end)^2), 1/delta_s(end)^2];
        end

    else %1<i<n
        switch order
            case 1
                nomin = [-delta_s(i), (delta_s(i) - delta_s(i-1)), delta_s(i-1)];
                denom = [delta_s(i-1)*(delta_s(i) + delta_s(i-1)), delta_s(i-1)*delta_s(i), delta_s(i)*(delta_s(i) + delta_s(i-1))];
                result(i, :) = circshift([nomin./denom, zeros(1, N-3)], i-2);
            case 2
                nomin = [2, -2, 2];
                denom = [delta_s(i-1)*(delta_s(i) + delta_s(i-1)), delta_s(i-1)*delta_s(i), delta_s(i)*(delta_s(i) + delta_s(i-1))];
                result(i, :) = circshift([nomin./denom, zeros(1, N-3)], i-2); 
        end
    end

end
