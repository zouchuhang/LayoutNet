function lines = compute_line_attributes(lines, vp)

%% determine if lines are above the horizon or below the horizon
[lines horizon] = is_line_above_horizon(lines, vp);
% disp_lines(rectimg, lines, [lines(:).above_horizon]);

%% determine if lines are left or right of the center vanishing point
[lines] = is_line_leftorright(lines, vp);
% disp_lines(rectimg, lines, [lines(:).leftorright]);

%% let lines carry their own ID
for i=1:length(lines), lines(i).id = i; end

%% compute 2D line equations for all lines
lines = compute_lineeq(lines);

%%
% lines = is_vertline_outside(lines);

