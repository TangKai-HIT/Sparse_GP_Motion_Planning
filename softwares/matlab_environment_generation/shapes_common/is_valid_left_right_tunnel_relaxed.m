%% 
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%%
function is_valid = is_valid_left_right_tunnel_relaxed( vertical_disp, bbox, gap_clearance)
%GET_WALL_WITH_UNIFORM_GAP Creates a wall with uniform gaps

is_valid = all(vertical_disp(1:(end-1)) >= bbox(3) + 0.51*gap_clearance ) & ...
            all(vertical_disp(1:(end-1)) <= bbox(4) - 0.51*gap_clearance ) & ...
            all(abs(diff(vertical_disp)) > 0.25*(bbox(4)-bbox(3)));
end