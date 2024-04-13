function res = quaterExpMap(quat)
%QUATEREXPMAP exponential map of a quaternion
%   此处显示详细说明
res = quat;

qv_norm = norm(quat(2:4));
res(1) = exp(quat(1)) * cos(qv_norm);
res(2:4) = exp(quat(1)) * sin(qv_norm) * quat(2:4) ./ qv_norm;
end

