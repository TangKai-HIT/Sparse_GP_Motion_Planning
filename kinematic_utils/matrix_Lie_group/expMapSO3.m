function SO3 = expMapSO3(omega, theta)
%EXPMAPSO3 

omegaNorm = norm(omega);
if omegaNorm>0
        omega = omega/omegaNorm;
end

if ~exist("theta","var")
    theta = omegaNorm;
end

omegaSkew = skewMatrix(omega);
SO3 = eye(3) + omegaSkew * sin(theta) + omegaSkew*omegaSkew*(1-cos(theta));