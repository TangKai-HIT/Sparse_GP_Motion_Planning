function [phi, axis] = unitQuat2rotVec(unit_quat)
%UNITQUAT2ROTVEC unit quaternion to rotation vector (through capitalized log map: S_3 -> R_3)
%   此处显示详细说明

qv_norm = norm(unit_quat(2:4));

phi = 2 * atan(qv_norm, unit_quat(1));

axis = unit_quat(2:4) ./ qv_norm;

end

