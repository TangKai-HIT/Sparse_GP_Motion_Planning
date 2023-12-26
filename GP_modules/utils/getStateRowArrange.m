function statesRow = getStateRowArrange(statesVec, stateDim)
%GETSTATEROWARRANGE turn states vector (N*stateDim) X 1 into numSupports X stateDim row array matrix
N = length(statesVec)/stateDim;
statesRow = reshape(statesVec, stateDim, N)';
end