function Phi = transMatMinAcc(t, s, dim)
%TRANSMATMINACC get 0-acceleration ode transition matrix
%   此处显示详细说明

I = eye(dim);
zero = zeros(dim,dim);
delta = t-s;
Phi = [I, delta*I;
      zero,  I];
end

