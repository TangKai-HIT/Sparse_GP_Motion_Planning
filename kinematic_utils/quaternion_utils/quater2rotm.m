function mat = quater2rotm(quat)
%QUATER2ROTM get equivalent matrix defining rotation by a unit quaternion 
%   此处显示详细说明

q_w = quat(1);
q_vec = zeros(3,1);
q_vec(1:3) = quat(2:4);

mat = (q_w^2 - q_vec' * q_vec) * eye(3) + 2 * (q_vec * q_vec') + 2 * q_w * skewMatrix(q_vec);
end

