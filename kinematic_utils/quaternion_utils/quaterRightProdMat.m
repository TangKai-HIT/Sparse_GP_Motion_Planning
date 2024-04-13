function mat = quaterRightProdMat(quat)
%QUATERRIGHTPRODMAT get right product matrix of a quaternion
%   此处显示详细说明

if size(quat, 2) ~=1
    quat = quat';
end

mat = quat(1) * eye(4) + [0, -quat(2:4)'; quat(2:4), -skewMatrix(quat(2:4))];
end

