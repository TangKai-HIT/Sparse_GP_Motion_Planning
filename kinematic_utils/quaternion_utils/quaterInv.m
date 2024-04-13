function quat_inv = quaterInv(quat)
%QUATERINV 此处显示有关此函数的摘要
%   此处显示详细说明

quat_inv = quaterConj(quat) ./ dot(quat, quat);
end

