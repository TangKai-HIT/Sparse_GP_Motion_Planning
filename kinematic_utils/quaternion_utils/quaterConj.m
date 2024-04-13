function quat_conj = quaterConj(quat)
%QUATERCONJ 此处显示有关此函数的摘要
%   此处显示详细说明

quat_conj = quat;
quat_conj(2:4) = -quat_conj(2:4);
end

