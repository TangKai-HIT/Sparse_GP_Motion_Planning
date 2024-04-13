function res = quaterCommut(quat1, quat2)
%QUATERPROD quaternion commutator
%   此处显示详细说明

res = zeros(size(quat1));

res(1) = 0;
res(2:4) = 2 .* cross(quat1(2:4), quat2(2:4));
end

