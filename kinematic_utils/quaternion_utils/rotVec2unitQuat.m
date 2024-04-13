function q = rotVec2unitQuat(phi, axis)
%ROTVEC2UNITQUAT rotation vector to unit quaternion (through capitalized exponential map: R_3 -> S_3)
%   此处显示详细说明

q = zeros(1, 4);
q(1) = cos(phi/2);
q(2:4) = sin(phi/2) .* axis;
end

