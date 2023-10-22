function Phi = transMatMinJerk(t, s, dim)
%TRANSMATMINJERK get 0-jerk ode transition matrix
%   此处显示详细说明

I = eye(dim);
zero = zeros(dim,dim);
delta = t-s; 
Phi = [I, delta.*I, 0.5*delta^2.*I;
          zero,  I,     delta.*I;
          zero, zero,       I];
end

