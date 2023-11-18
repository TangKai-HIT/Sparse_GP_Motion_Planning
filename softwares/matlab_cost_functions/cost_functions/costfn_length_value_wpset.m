%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function value = costfn_length_value_wpset( wpset )
%COSTFN_LENGTH_VALUE_WPSET Get length of wpset trajectory
%   wpset: NxM trajectory
%   value: length of trajectory

dxi = diff(wpset,1,1);
value = sum(sqrt(sum(dxi.^2,2)));

end

