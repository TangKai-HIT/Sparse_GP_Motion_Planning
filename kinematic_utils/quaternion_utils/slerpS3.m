function quat_interp = slerpS3(quat0, quat1, t_squence)
%SLERPS3 spherical linear interpolation on S3
%   inputs:
%       quat0, quat1: unit quaternion, elements on S3
%       t_squence: 1 X N or N x 1, ranging in [0, 1]
%   outputs:
%       quat_interp: N X 4

cos_theta = dot(quat0, quat1);

% check if angle =< pi/2 (rotation angle <= pi), ensure shortest path interpolation
if cos_theta>=0
    del_theta = acos(cos_theta);
else
    quat1 = -quat1;
    del_theta = acos(dot(quat0, quat1));
end

if size(quat0, 1) ~=1
    quat0 = quat0';
end

if size(quat1, 1) ~=1
    quat1 = quat1';
end

if size(t_squence, 2) ~=1
    t_squence = t_squence';
end

quat_interp = sin((1-t_squence) * del_theta)./sin(del_theta) * quat0 + sin(t_squence * del_theta)./sin(del_theta) * quat1;

end

