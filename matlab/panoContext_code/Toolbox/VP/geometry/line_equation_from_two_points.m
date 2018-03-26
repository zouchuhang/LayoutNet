function lineeq = line_equation_from_two_points(p1, p2)
% function lineeq = line_equation_from_two_points(p1, p2)
% returns line parameters(ax+by+c=0) that passes through p1 and p2
% p1, p2: [2x1] (x,y)
% lineeq: [3x1] (a,b,c)

lineeq = [p2(2)-p1(2) p1(1)-p2(1) -p1(1)*p2(2)+p2(1)*p1(2)]';
lineeq = lineeq / norm(lineeq(1:2));

% normalize sign to be consistent... (0,0) should have a positive distance
distzero = [0 0 1] * lineeq;
if abs(distzero) > 1e-10
    lineeq = lineeq * sign(distzero);
end

