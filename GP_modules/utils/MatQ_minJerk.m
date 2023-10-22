function Q = MatQ_minJerk(a, b, Q_c)
%MATQ_MINACC Q matrix for 0-jerk ode, Q_a_b = intergrate Phi(b,s)*F(s)*Q_c*F(s)'*Phi(b,s)'*ds from a to b,
%   where F(s) = [0; I]
%   此处显示详细说明

delta = b-a;
Q = [1/2*delta^5*Q_c,   1/8*delta^4*Q_c,    1/6*delta^3*Q_c;
        1/8*delta^4*Q_c,    1/3*delta^3*Q_c,    1/2*delta^2*Q_c;
        1/6*delta^3*Q_c,    1/2*delta^2*Q_c,    delta*Q_c];
end

