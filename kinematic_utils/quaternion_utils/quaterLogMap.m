function res = quaterLogMap(quat)
%QUATERLOGMAP logarithm of a quaternion
%   此处显示详细说明
res = quat;

qv_norm = norm(quat(2:4));
theta = atan(qv_norm, quat(1));
u_vec = quat(2:4) / qv_norm;
res(1) = log(qv_norm);
res(2:4) = u_vec * theta;
end

