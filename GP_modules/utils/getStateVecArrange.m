function statesVec = getStateVecArrange(statesRow, stateDim)
%GETSTATEVECARRANGE turn N X stateDim row array matrix into states vector (N*stateDim) X 1
statesVec = reshape(statesRow', size(statesRow, 1)*stateDim, 1);
end