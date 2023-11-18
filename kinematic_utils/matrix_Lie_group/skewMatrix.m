function skew = skewMatrix(vec)

N = length(vec);

switch N
    case 1
        skew = [0  -vec;
                      vec  0];
    case 3
        skew = [0      -vec(3)  vec(2);
                     vec(3)     0       -vec(1);
                     -vec(2)    vec(1)       0];
    otherwise
        skew = [];
        disp("R1-scalar or R3-vector");
end