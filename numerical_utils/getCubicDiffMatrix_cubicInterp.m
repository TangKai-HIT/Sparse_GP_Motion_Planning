function result = getCubicDiffMatrix_cubicInterp(points)
%GETCUBICDIFFMATRIX_CUBICINTERP get 3rd differential matrix of unevenly sampled points using cubic newton-form polynomial interpolation
%                                                           derived from "table of divided difference"  
%   Inputs:
%       points: N X 1 or 1 X N, sampled points
%   Output:
%       result: N X N

N = length(points);
delta_s = diff(points, 1);
result = zeros(N, N);

diff2_mat = getDiffMatrix_quadInterp(points, 2); %use for boundary condition

for i = 1:N
    if i==1
        result(i, :) = (diff2_mat(2,:) - diff2_mat(1,:)) ./ delta_s(1);

    elseif i==N-1
        result(i, :) = (diff2_mat(N-1,:) - diff2_mat(N-2,:)) ./ delta_s(end-1);
    
    elseif i==N
        result(i, :) = (diff2_mat(N,:) - diff2_mat(N-1,:)) ./ delta_s(end);

    else %1<i<n-1
        ds_0 = delta_s(i-1);
        ds_1 = delta_s(i);
        ds_2 = delta_s(i+1);
        
        cur_term =  [1/(ds_2*(ds_1 + ds_2)*(ds_0 + ds_1 + ds_2)), -1/(ds_1*ds_2*(ds_0 + ds_1)), 1/(ds_0*ds_1*(ds_1 + ds_2)), -1/(ds_0*(ds_0 + ds_1)*(ds_0 + ds_1 + ds_2)), zeros(1, N-4)];
        result(i, :) = circshift(cur_term, i-2); 
        
    end

end
