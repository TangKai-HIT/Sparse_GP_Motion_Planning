function plotViaPtSpline1D(t_sample, q, dq, ddq, dddq, t_wp, q_wp)
%PLOTVIAPTSPLINE1D 此处显示有关此函数的摘要
%   此处显示详细说明

figure('Position',[500,100, 1000,800]);
subplot(2,2,1)
plot(t_sample, q, '-r', 'LineWidth', 0.8); hold on;
if exist("t_wp", "var") && exist("q_wp", "var")
    scatter(t_wp, q_wp, 'black', 'filled' ,'o', SizeData=20);
end
title('pos');

subplot(2,2,2)
plot(t_sample, dq, '-r', 'LineWidth', 0.8); 
title('vel');

subplot(2,2,3)
plot(t_sample, ddq, '-r', 'LineWidth', 0.8); 
title('accel');

subplot(2,2,4)
plot(t_sample, dddq, '-r', 'LineWidth', 0.8); 
title('jerk');
end

