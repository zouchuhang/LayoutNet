function [lines horizon] = is_line_leftorright(lines, vp)
%
% left or right of vp{3}.. vp{3} is the vp close to image center
% adds lines(i).leftorright:  1 if right, -1 if left, 0 if neither
% horizon: [3x1] (a,b,c) parameters of the equation of the horizon in 2D where ax+by+c=0

%% centerline is the line connecting vp{1} and vp{3}
centerline = line_equation_from_two_points(vp{1}, vp{3});

if distance_of_point_to_line(centerline, [10000000 0]) < 0 % if a really right point is negative
	centerline = -centerline; % flip sign of horizon
end

%% test all lines
for i = 1:length(lines)
	l1 = sign(distance_of_point_to_line(centerline, lines(i).point1));
	l2 = sign(distance_of_point_to_line(centerline, lines(i).point2));
	if l1==1 && l2==1
		lines(i).leftorright = 1;
	elseif l1==-1 && l2==-1
		lines(i).leftorright = -1;
	else
		lines(i).leftorright = 0;
	end
end

