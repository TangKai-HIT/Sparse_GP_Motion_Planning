function [t_l, t_u, l_id, u_id] = findInterval(t, t_sq)
%FINDINTERVAL: find interval of t in a ascending sequence
%   Inputs:
%       t_sq: size>=2
%   Output:
%       t_l, t_u: lower & upper bound of interval, if exceed (t_u=[]), if beyond (t_l=[])
%       l_id, u_id: index of interval's lower & upper bound

N = length(t_sq);

l = 0; u = N+1;
m = floor((l+u)/2);

while t<t_sq(m) || t>t_sq(m+1)
    if t<t_sq(m)
        u = m;
    elseif t>t_sq(m+1)
        l = m+1;
    end

    m = floor((l+u)/2);
    if m==0 || m==N
        break;
    end
end

if m<1 %t<a, t_sq in [a,b]
    t_l = [];
    t_u = t_sq(1);
    l_id = []; 
    u_id = 1;
elseif m>=N
    if t>t_sq(m) %t > b
        t_l = t_sq(N);
        t_u = [];
        l_id = N; 
        u_id = [];
    else %t==b
        t_l = t_sq(m-1);
        t_u = t_sq(m);
        l_id = m-1; 
        u_id = m;
    end
else %[l, u)
    if t==t_sq(m+1) && m+1<N
        m=m+1;
    end
    t_l = t_sq(m);
    t_u = t_sq(m+1);
    l_id = m; 
    u_id = m+1;
end

end