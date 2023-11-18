function  jacobRInv = jacobRInvSO3(theta)
%JACOBRSO3  inverse of right Jacobian of SO3
%   theta: vedge operation from a se3

thetaNorm = norm(theta);
thetaSkew = skewMatrix(theta);

jacobRInv = eye(3) + 0.5*thetaSkew + (1/thetaNorm^2 - (1+cos(thetaNorm))/(2*thetaNorm*sin(thetaNorm))) * (thetaSkew * thetaSkew);
end

