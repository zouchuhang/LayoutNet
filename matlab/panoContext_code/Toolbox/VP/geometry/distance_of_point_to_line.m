function d = distance_of_point_to_line(line, x)
% function d = distance_of_point_to_line(line, x)
% 
% returns d = ax+by+c
% line: [3x1] (a,b,c)
% x: [2x1] (x,y)

xh = [x(:); 1];
d = line(:)' * xh;
