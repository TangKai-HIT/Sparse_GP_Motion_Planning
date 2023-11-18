%%
% Copyright (c) 2016 Carnegie Mellon University, Sanjiban Choudhury <sanjibac@andrew.cmu.edu>
%
% For License information please see the LICENSE file in the root directory.
%
%% EXAMPLE_TUNNEL
% Make a left to right tunnel
clc;
clear;
close all;

%% Create two rectangles
scenario = 3;
switch scenario
    case 1
        bbox = [0 1 0 1]; %unit bounding box
        vertical_disp = [0.4 0.7 0.2 0.5 1.0];
        wall_width = 0.02;
        gap_clearance = 0.1;
        rectangle_array = get_left_right_tunnel( bbox, vertical_disp, wall_width, gap_clearance);
    case 2
        rng(5);
        bbox = [0 1 0 1]; %unit bounding box
        wall_width = 0.02;
        gap_clearance = 0.1;
        while true
            vertical_disp = transpose(cumsum(randfixedsum(5,1,1.0,-0.5,0.5)));
            if(is_valid_left_right_tunnel( vertical_disp, bbox, gap_clearance))
                break;
            end
        end
        rectangle_array = get_left_right_tunnel( bbox, vertical_disp, wall_width, gap_clearance);
%     case 3
%         bbox = [0 1000 -100 100]; %unit bounding box
%         vertical_disp = [0 -30 20 -30 20];
%         wall_width = 1;
%         gap_clearance = 50;
%         rectangle_array = get_left_right_tunnel( bbox, vertical_disp, wall_width, gap_clearance);
    case 3
        bbox = [0 1000 0 1000]; %unit bounding box
        wall_width = 10;
        gap_clearance = 50.0;
        while true
            vertical_disp = transpose(cumsum( 1000*(rand(20,1) - 0.5) ));
            %vertical_disp = transpose(cumsum(randfixedsum(10,1, 1.0,-0.5,0.5)));
            if(is_valid_left_right_tunnel_relaxed( vertical_disp, bbox, gap_clearance))
                break;
            end
        end

        rectangle_array = get_left_right_tunnel( bbox, vertical_disp, wall_width, gap_clearance);
        bbox = [-200 1200 -200 1200];
        
%         x = -gap_clearance*0.5 + gap_clearance*(1:length(vertical_disp));
%         coord = [x' vertical_disp' [x(2:end) 1200]' vertical_disp'];
%         coord = reshape(coord', 2 , [])';
%         coord = [gap_clearance*0.5 0; coord];
end

%% Visualize
figure;
visualize_shapes(rectangle_array);
axis equal
axis(bbox)

