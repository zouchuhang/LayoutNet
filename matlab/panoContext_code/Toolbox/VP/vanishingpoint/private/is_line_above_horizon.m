function [lines horizon] = is_line_above_horizon(lines, vp)
%
% adds lines(i).above_horizon:  1 if above, -1 if below, 0 if neither
% horizon: [3x1] (a,b,c) parameters of the equation of the horizon in 2D where ax+by+c=0

%% horizon is the line connecting vp{2} and vp{3}
horizon = line_equation_from_two_points(vp{2}, vp{3});

if distance_of_point_to_line(horizon, [0 -10000000]) < 0 % if a really high point is not above horizon
	horizon = -horizon; % flip sign of horizon
end

%% test all lines
for i = 1:length(lines)
	l1 = sign(distance_of_point_to_line(horizon, lines(i).point1));
	l2 = sign(distance_of_point_to_line(horizon, lines(i).point2));
	if l1==1 && l2==1
		lines(i).above_horizon = 1;
	elseif l1==-1 && l2==-1
		lines(i).above_horizon = -1;
	else
		lines(i).above_horizon = 0;
	end
end

