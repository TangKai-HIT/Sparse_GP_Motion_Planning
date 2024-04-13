function res = quaterRotVec(vec, quat)
%QUATERROTVEC rotate a vector by a unit quaternion
%   此处显示详细说明

vec_quat = zeros(size(quat));
vec_quat(2:4) = vec;

quat_conj = quaterConj(quat);
res = quaterProd(quaterProd(quat, vec_quat), quat_conj);
res = res(2:4);
end

