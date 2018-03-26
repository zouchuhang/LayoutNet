function [regionid boundary] = getregionid(x,y,vp)

% input: x [n x 1], y [n x 1], vp
% output: regionid[n x 1] {1, 2, 3, 4}
%  1 | 2
% ---+---
%  3 | 4

%%
p = [x(:) y(:) ones(length(x), 1)];


%% horizon is the line connecting vp{2} and vp{3}
horizon = line_equation_from_two_points(vp{2}, vp{3});
if distance_of_point_to_line(horizon, [0 -10000000]) < 0 % if a really high point is not above horizon
	horizon = -horizon; % flip sign of horizon
end

%% centerline is the line connecting vp{1} and vp{3}
centerline = line_equation_from_two_points(vp{1}, vp{3});
if distance_of_point_to_line(centerline, [10000000 0]) < 0 % if a really right point is negative
	centerline = -centerline; % flip sign of horizon
end

d1 = horizon(:)' * p';
ud = (d1 < 0);   % 0 if up, 1 if down
    
d2 = centerline(:)' * p';
lr = (d2 > 0);    % 0 if left, 1 if right

regionid = ud*2 + lr + 1;
%  1 | 2
% ---+---
%  3 | 4

%%
if nargout>1
    BOUNDARY_DIST_THRES = 15;
    boundary = 0;
    
    if abs(d1)/norm(horizon(1:2)) < BOUNDARY_DIST_THRES || ...
       abs(d2)/norm(centerline(1:2)) < BOUNDARY_DIST_THRES
        boundary = 1;
    end
end
