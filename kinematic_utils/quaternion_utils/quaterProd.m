function res = quaterProd(quat1, quat2)
%QUATERPROD quaternion product
%   此处显示详细说明

res = zeros(size(quat1));

res(1) = quat1(1)*quat2(1) - dot(quat1(2:4), quat2(2:4));
res(2:4) = quat1(1)*quat2(2:4) + quat2(1)*quat1(2:4) + cross(quat1(2:4), quat2(2:4));
end

