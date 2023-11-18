function  jacobR = jacobRSO3(theta)
%JACOBRSO3 right Jacobian of SO3
%   theta: vedge operation from a se3

thetaNorm = norm(theta);
thetaSkew = skewMatrix(theta);

jacobR = eye(3) - (1 - cos(thetaNorm))*thetaSkew/(thetaNorm^2) + (thetaNorm - sin(thetaNorm))* (thetaSkew * thetaSkew)/(thetaNorm^3);
end

