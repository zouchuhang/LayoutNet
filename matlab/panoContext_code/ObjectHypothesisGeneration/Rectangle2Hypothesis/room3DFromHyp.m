function [ objects] = room3DFromHyp( points, floor_height )
%ROOM3DFROMHYP Popup room layout hypothesis to 3D
%   points: room corner in sequence in room layout
%   floor_height: assumed floor height, e.g. -160

% uv = xyz2uvN(points);
% high_height = floor_height * tan(uv(1:4,2))./tan(uv(5:8,2));
% ceil_height = mean(high_height);
% 
% point3D = zeros(8,3);
% for i = 5:8
%     point3D(i,:) = LineFaceIntersection([0 0 floor_height], [0 0 1], [0 0 0], points(i,:));
% end
COUNTER_CLOCKWISE = true;
FITTINGRULES;

[pt I] = sortXYZ(points(1:4,:));
pd = points(5:8,:); pd = pd(I,:);
sort_points = [pd;pt];
ann_points = sort_points(POINTMAPPING{3},:);

objects = struct('out_points_w',[],'name',[],'type',[],'points',[],'x_w',[]);
[ out_point, x, fval, succ ] = initial_cuboid_3D( ann_points, 3, floor_height, 3, true );

objects.points = ann_points;
objects.type = 3;
objects.name = 'room';
objects.out_points_w = out_point;
objects.x_w = x;

end

