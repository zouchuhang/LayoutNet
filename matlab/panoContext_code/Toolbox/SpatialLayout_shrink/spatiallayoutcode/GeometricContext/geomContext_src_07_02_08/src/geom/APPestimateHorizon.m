function y_h = APPestimateHorizon(lines)
% Estimates horizon position of image
% lines should be normalized according to size of image
% need to multiply output by mean(imsize(1:2)) to get pixel value
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005

max_angle = pi/8;

% lines(num_lines, [x1 x2 y1 y2 angle r])
% points(num_points, [y x])


a = lines(:, 5);
x1 = lines(:, 1);
y1 = lines(:, 3);
x2 = lines(:, 2);
y2 = lines(:, 4);
r = lines(:, 6);
cosa = cos(a);
sina= sin(a);
tana = tan(a);

max_dist = 1E2;

num_lines = size(lines, 1);
if num_lines <= 1
    points = zeros(0, 2);
    return;
end

num_pairs = num_lines*(num_lines-1)/2;

adist = pdist(a, 'cityblock');
%adist = mod(adist, pi/2) < max_angle;
adist = adist < max_angle;
adistM = squareform(adist);

points = zeros(num_pairs, 2);
count = 0;
for i = 1:num_lines
    for j = i+1:num_lines
        if adistM(i, j) & (r(i)~=r(j) | a(i)~=a(j))
            count = count+1;            
            if a(i) == a(j) | a(i) == -a(j) % parallel but not colinear
                points(count, 1:2) = ([x1(i) y1(i)]+[x1(j) y1(j)])/2 + [sina(i) cosa(i)]*max_dist;
            else
                m1 = tana(i);
                m2 = tana(j);
                points(count, 1) = (-y1(i)*m2+m1*y1(j)-m1*m2*x1(j)+m1*x1(i)*m2)/(m1-m2);
                points(count, 2) = (y1(j)-y1(i)+m1*x1(i)-m2*x1(j))/(m1-m2);                
                if abs(points(count, 1))==Inf | abs(points(count, 2))==Inf
                    points(count, 1:2) = ([x1(i) y1(i)]+[x1(j) y1(j)])/2 + [sina(i) cosa(i)]*max_dist;
                end
                    
            end
        end
    end
end

points = points(1:count, :);

%ind = find(abs(points(:, 2))<0.02);% & abs(points(:, 1)) < 0.75);
%tpoints = points(ind, :);

if size(points, 1) == 0
    y_h = 0.5;
else
    y_h = 0.5+fminbnd(@sum_distance, -2, 2, [], points(:, 1));
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = sum_distance(x, pts)
d = sum(sqrt(abs(x-pts)));
