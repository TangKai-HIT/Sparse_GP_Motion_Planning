function vec = skew2Vec(skew)

N = size(skew, 1);

switch N
    case 2
        vec = skew(2,1);
    case 3
        vec = [skew(3,2); -skew(3,1); skew(2,1)];
    otherwise
        vec = [];
end