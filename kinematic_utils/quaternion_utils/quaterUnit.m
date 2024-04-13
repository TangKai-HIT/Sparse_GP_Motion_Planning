function quat = quaterUnit(theta, u_vec)
%QUATERUNIT set unit quaternion
%   此处显示详细说明
quat = [cos(theta), u_vec*sin(theta)];
end

