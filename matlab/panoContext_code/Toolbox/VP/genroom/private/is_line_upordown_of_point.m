function ud = is_line_upordown_of_point(lines, point)
% ud:  1 if up, -1 if down, 0 if neither

p1 = cat(1, lines.point1);
p2 = cat(1, lines.point2);

ud1 = p1(2) < point(2);
ud2 = p2(2) < point(2);

ud = 1 * (ud1==1 & ud2==1) + (-1) * (ud1==0 & ud2==0);
