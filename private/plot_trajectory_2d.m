function fig = plot_trajectory_2d(results)
%PLOT_TRAJECTORY_2D 显示全局二维位移轨迹，标记起点、终点及模板更新帧。
% 使用原结果的X/Y分量，单位nm，不交换坐标、不平滑轨迹。
xy = results.traj_global;
fig = figure('Name', '2D displacement trajectory', 'Color', 'w');
ax = axes('Parent', fig);
hold(ax, 'on');
plot(ax, xy(:,1), xy(:,2), '.-', 'Color', [0.15 0.4 0.7], ...
    'LineWidth', 1.2, 'MarkerSize', 8, 'DisplayName', 'Trajectory');
plot(ax, xy(1,1), xy(1,2), 'o', 'MarkerSize', 8, ...
    'MarkerEdgeColor', [0.1 0.55 0.2], 'MarkerFaceColor', [0.1 0.55 0.2], ...
    'DisplayName', 'Start');
plot(ax, xy(end,1), xy(end,2), 's', 'MarkerSize', 8, ...
    'MarkerEdgeColor', [0.8 0.15 0.15], 'MarkerFaceColor', [0.8 0.15 0.15], ...
    'DisplayName', 'End');
if ~isempty(results.updateFrames)
    idx = results.updateFrames;
    plot(ax, xy(idx,1), xy(idx,2), 'd', 'MarkerSize', 7, ...
        'MarkerEdgeColor', [0.85 0.5 0.05], 'MarkerFaceColor', [1 0.75 0.2], ...
        'DisplayName', 'Template update');
end
xlabel(ax, 'X displacement (nm)');
ylabel(ax, 'Y displacement (nm)');
title(ax, 'Global 2D displacement trajectory');
axis(ax, 'equal');
% 留出边距，使起点和终点标记不会被坐标轴边界裁切。
lo = min(xy, [], 1);
hi = max(xy, [], 1);
pad = 0.05 * max([hi - lo, 1]);
xlim(ax, [lo(1) - pad, hi(1) + pad]);
ylim(ax, [lo(2) - pad, hi(2) + pad]);
grid(ax, 'on');
box(ax, 'on');
legend(ax, 'show', 'Location', 'northwest');
hold(ax, 'off');
end
