function [omega, theta] = logMapSO3(SO3)

if trace((SO3-eye(3))'*(SO3-eye(3)))==0
    theta = 0;
    omega = [1;0;0];

elseif trace(SO3) == -1
    theta = pi;
    omega = 1/sqrt(2*(1+SO3(3,3))) *  [SO3(1,3); SO3(2,3); 1+SO3(3,3)];

else
    theta = acos(0.5*(trace(SO3) - 1));
    omegaSkew = (SO3 - SO3') / (2*sin(theta));
    omega = skew2Vec(omegaSkew);

end

if nargout==1
    omega = omega * theta;
end