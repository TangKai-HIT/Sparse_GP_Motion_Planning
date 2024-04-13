function quat_interp = slerpS3LogMap(quat0, quat1, t_squence)
%SLERPS3LOGMAP spherical linear interpolation on S3 using Log map (low efficiency, use slerpS3 instead !)
%   inputs:
%       t_squence: 1 X N or N x 1, ranging in [0, 1]

del_quat = quaterProd(quaterConj(quat0), quat1);
[phi, axis] = unitQuat2rotVec(del_quat);

leftMat0 = quaterLeftProdMat(quat0);

if size(axis, 1) ~=1
    axis = axis';
end

if size(t_squence, 2) ~=1
    t_squence = t_squence';
end

quat_interp = leftMat0 * [cos(phi/2 .* t_squence),  sin(phi/2 .* t_squence) * axis];

end

