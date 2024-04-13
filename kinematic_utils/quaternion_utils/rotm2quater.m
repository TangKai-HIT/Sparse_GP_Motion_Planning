function quat = rotm2quater(rotm)
%ROTM2QUATER SO3 to unit quaternion 
%   此处显示详细说明

q_w = 0.5 * sqrt(1 + trace(rotm));
q_v = 1/(4*q_w) * [rotm(3,2)-rotm(2,3),  rotm(1,3)-rotm(3,1),  rotm(2,1)-rotm(1,2)];

quat = [q_w, q_v];
end

