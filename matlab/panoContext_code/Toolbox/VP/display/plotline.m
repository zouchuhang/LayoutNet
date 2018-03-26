function h = plotline(p1, p2, color, linewidth)
if ~exist('color', 'var')
	color = [0 1 0];
end
if ~exist('linewidth','var')
    linewidth = 1;
end

h = plot([p1(1) p2(1)], [p1(2) p2(2)], 'Color', color, 'LineWidth', linewidth);
