function [ objects ] = rectangleCombo2Hypothesis( xyzBoxes, views, room_points )
%RECTANGLECOMBO2HYPOTHESIS Summary of this function goes here
%   Detailed explanation goes here
COUNTER_CLOCKWISE = true;
FITTINGRULES;
uni_room_points = zeros(8,3);
uni_room_points(POINTMAPPING{3},1) = room_points(:,1);
uni_room_points(POINTMAPPING{3},2) = room_points(:,2);
uni_room_points(POINTMAPPING{3},3) = room_points(:,3);

I = find(uni_room_points(1:4,1)<0 & uni_room_points(1:4,2)<0);
ID = rem([I I+1 I+2 I+3]-1,4) + 1;
uni_room_points = uni_room_points([ID ID+4],:);
wall_points = zeros(4,3,6);
wall_points(:,:,1) = uni_room_points([1 4 8 5],:);
wall_points(:,:,2) = uni_room_points([3 2 6 7],:);
wall_points(:,:,3) = uni_room_points([2 1 5 6],:);
wall_points(:,:,4) = uni_room_points([4 3 7 8],:);
wall_points(:,:,5) = uni_room_points([1 2 3 4],:);
wall_points(:,:,6) = uni_room_points([8 7 6 5],:);

LRDU_SURF = zeros(6,4);
LRDU_SURF(1,:) = [3 4 5 6];
LRDU_SURF(2,:) = [4 3 5 6];
LRDU_SURF(3,:) = [2 1 5 6];
LRDU_SURF(4,:) = [1 2 5 6];
LRDU_SURF(5,:) = [1 2 4 3];
LRDU_SURF(6,:) = [1 2 3 4];

normals = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
normal_ids = [1 1 2 2 3 3];
main_angle = [pi/2 -pi/2 -pi 0 0];

objects = repmat(struct('out_points_w',[],'name',[],'type',[],'points',[],'x_w',[]),0,1);
object = struct('out_points_w',[],'name','other','type',[],'points',[],'x_w',zeros(7,1));
xyzBoxes = reshape(xyzBoxes', 3, 8)';

I = views(1)==[1;4;2;3] & views(2)==[4;2;3;1];
if any(I)
    t = views(2);
    views(2) = views(1);
    views(1) = t;
    t = xyzBoxes(1:4,:);
    xyzBoxes(1:4,:) = xyzBoxes(5:8,:);
    xyzBoxes(5:8,:) = t;
end

xyzBoxes_norm = xyzBoxes./repmat(sqrt(sum(xyzBoxes.^2,2)),1,3);
uvBoxes = xyz2uvN(xyzBoxes_norm);

%% on floor?
[inside,~,~] = insideCone( wall_points(:,:,5), xyzBoxes, 0 );
if any(inside)
    lhs_gnd_pt1 = LineFaceIntersection(uni_room_points(1,:), [0 0 1], [0 0 0], xyzBoxes(4,:));
    lhs_gnd_pt2 = LineFaceIntersection(uni_room_points(1,:), [0 0 1], [0 0 0], xyzBoxes(3,:));
    rhs_gnd_pt1 = LineFaceIntersection(uni_room_points(1,:), [0 0 1], [0 0 0], xyzBoxes(8,:));
    rhs_gnd_pt2 = LineFaceIntersection(uni_room_points(1,:), [0 0 1], [0 0 0], xyzBoxes(7,:));
    
    point2 = (lhs_gnd_pt2 + rhs_gnd_pt1)/2;
    point3 = point2;    point6 = point2;
    d3 = setdiff([1 2], normal_ids(views(1)));
    d6 = setdiff([1 2], normal_ids(views(2)));
    point3(d3) = lhs_gnd_pt1(d3);
    point6(d6) = rhs_gnd_pt2(d6);
    
    camera_height = uni_room_points(1,3);
    h1 = camera_height/tan(uvBoxes(4,2))*tan(uvBoxes(1,2));
    h2 = camera_height/tan(uvBoxes(3,2))*tan(uvBoxes(2,2));
    h3 = camera_height/tan(uvBoxes(8,2))*tan(uvBoxes(5,2));
    h4 = camera_height/tan(uvBoxes(7,2))*tan(uvBoxes(6,2));
    ceil_z = (h1 + h2 + h3 + h4)/4;
    point1 = point2;    point4 = point3;    point5 = point6;
    point1(3) = ceil_z;
    point4(3) = ceil_z;
    point5(3) = ceil_z;
    
    object.out_points_w = [point1; point2; point3; point4; point5; point6];
%     object.type = 9;
%     object.points = object.out_points_w ./ repmat(sqrt(sum(object.out_points_w.^2,2)),1,3);
%     object.x_w(1:3) = point3;
%     object.x_w(4) = norm(point2-point3);
%     object.x_w(5) = norm(point6-point2);
%     object.x_w(6) = norm(point4-point3);
%     object.x_w(7) = main_angle(views(1));
    objects(end+1,1) = object;
end

%% on lhs wall of the lhs views
[inside,~,~] = insideCone( wall_points(:,:,LRDU_SURF(views(1),1)), xyzBoxes(1:4,:), 0 );
if any(inside)
    LID = LRDU_SURF(views(1),1);
    wall_point = wall_points(:,:,LID);
    
    point4 = LineFaceIntersection(wall_point(1,:), normals(LID,:), [0 0 0], xyzBoxes(1,:));
    point3 = LineFaceIntersection(wall_point(1,:), normals(LID,:), [0 0 0], xyzBoxes(4,:));
    
    point1_1 = LineFaceIntersection(point4, normals(views(1),:), [0 0 0], xyzBoxes(5,:));
    point2_1 = LineFaceIntersection(point3, normals(views(1),:), [0 0 0], xyzBoxes(8,:));
%     point5_1 = LineFaceIntersection(point1_1, normals(views(2),:), [0 0 0], xyzBoxes(6,:));
%     point6_1 = LineFaceIntersection(point2_1, normals(views(2),:), [0 0 0], xyzBoxes(7,:));
    
    point1_2 = LineFaceIntersection(point4, normals(views(1),:), [0 0 0], xyzBoxes(2,:));
    point2_2 = LineFaceIntersection(point3, normals(views(1),:), [0 0 0], xyzBoxes(3,:));  
    
    id = normal_ids(LID);
    point1 = point4;    point2 = point3;
    val = (point1_1(id) + point1_2(id) + point2_1(id) + point2_2(id))/4;
    point1(id) = val; %(point1_1(id) + point1_2(id))/2;
    point2(id) = val; %(point2_1(id) + point2_2(id))/2;
    
    point5_t = LineFaceIntersection(point1, normals(views(2),:), [0 0 0], xyzBoxes(6,:));
    point6_t = LineFaceIntersection(point2, normals(views(2),:), [0 0 0], xyzBoxes(7,:));
    id = setdiff([1 2],id);
    point5 = point1;
    point6 = point2;
    point5(id) = point5_t(id);
    point6(id) = point6_t(id);
    
%     low_h = (point3(3)+point6_1(3))/2;
%     high_h = (point4(3)+point5_1(3))/2; 
%     left_xy = ( point3(1:2) + point4(1:2) )/2;
%     middle_xy = ( point1_1(1:2) + point2_1(1:2))/2;
%     right_xy = (point5_1(1:2) + point6_1(1:2))/2;
%     
%     object.out_points_w = [middle_xy high_h; ...
%                            middle_xy low_h; ...
%                            left_xy low_h; ...
%                            left_xy high_h; ...
%                            right_xy high_h; ...
%                            right_xy low_h];
%     objects(end+1,1) = object;
    
%     low_h = (point3(3)+point6_2(3))/2;
%     high_h = (point4(3)+point5_2(3))/2; 
%     left_xy = ( point3(1:2) + point4(1:2) )/2;
%     middle_xy = ( point1_2(1:2) + point2_2(1:2))/2;
%     right_xy = (point5_2(1:2) + point6_2(1:2))/2;
%     
%     object.out_points_w = [middle_xy high_h; ...
%                            middle_xy low_h; ...
%                            left_xy low_h; ...
%                            left_xy high_h; ...
%                            right_xy high_h; ...
%                            right_xy low_h];
    object.out_points_w = [point1; point2; point3; point4; point5; point6];
    objects(end+1,1) = object; 
     
end

%% on rhs wall of the rhs views
[inside,~,~] = insideCone( wall_points(:,:,LRDU_SURF(views(2),2)), xyzBoxes(5:8,:), 0 );
if any(inside)
    RID = LRDU_SURF(views(2),2);
    wall_point = wall_points(:,:,RID);
    
    point5 = LineFaceIntersection(wall_point(1,:), normals(RID,:), [0 0 0], xyzBoxes(6,:));
    point6 = LineFaceIntersection(wall_point(1,:), normals(RID,:), [0 0 0], xyzBoxes(7,:));
    
    point1_1 = LineFaceIntersection(point5, normals(views(2),:), [0 0 0], xyzBoxes(2,:));
    point2_1 = LineFaceIntersection(point6, normals(views(2),:), [0 0 0], xyzBoxes(3,:));
%     point4_1 = LineFaceIntersection(point1_1, normals(views(1),:), [0 0 0], xyzBoxes(1,:));
%     point3_1 = LineFaceIntersection(point2_1, normals(views(1),:), [0 0 0], xyzBoxes(4,:));
    
    point1_2 = LineFaceIntersection(point5, normals(views(2),:), [0 0 0], xyzBoxes(5,:));
    point2_2 = LineFaceIntersection(point6, normals(views(2),:), [0 0 0], xyzBoxes(8,:));
    
    id = normal_ids(RID);
    point1 = point5;    point2 = point6;
    val = (point1_1(id) + point1_2(id) + point2_1(id) + point2_2(id))/4;
    point1(id) = val; %(point1_1(id) + point1_2(id))/2;
    point2(id) = val; %(point2_1(id) + point2_2(id))/2;
    
    point4_t = LineFaceIntersection(point1, normals(views(1),:), [0 0 0], xyzBoxes(1,:));
    point3_t = LineFaceIntersection(point2, normals(views(1),:), [0 0 0], xyzBoxes(4,:));
    id = setdiff([1 2],id);
    point4 = point1;
    point3 = point2;
    point4(id) = point4_t(id);
    point3(id) = point3_t(id);
    
    
%     point4_2 = LineFaceIntersection(point1_2, normals(views(1),:), [0 0 0], xyzBoxes(1,:));
%     point3_2 = LineFaceIntersection(point2_2, normals(views(1),:), [0 0 0], xyzBoxes(4,:));
    
%     low_h = (point3_1(3)+point6(3))/2;
%     high_h = (point4_1(3)+point5(3))/2; 
%     left_xy = ( point3_1(1:2) + point4_1(1:2) )/2;
%     middle_xy = ( point1_1(1:2) + point2_1(1:2))/2;
%     right_xy = (point5(1:2) + point6(1:2))/2;
%     
%     object.out_points_w = [middle_xy high_h; ...
%                            middle_xy low_h; ...
%                            left_xy low_h; ...
%                            left_xy high_h; ...
%                            right_xy high_h; ...
%                            right_xy low_h];
%     objects(end+1,1) = object;
    
%     low_h = (point3_2(3)+point6(3))/2;
%     high_h = (point4_2(3)+point5(3))/2; 
%     left_xy = ( point3_2(1:2) + point4_2(1:2) )/2;
%     middle_xy = ( point1_2(1:2) + point2_2(1:2))/2;
%     right_xy = (point5(1:2) + point6(1:2))/2;
%     
%     object.out_points_w = [middle_xy high_h; ...
%                            middle_xy low_h; ...
%                            left_xy low_h; ...
%                            left_xy high_h; ...
%                            right_xy high_h; ...
%                            right_xy low_h];
    object.out_points_w = [point1; point2; point3; point4; point5; point6];
    objects(end+1,1) = object;     
end

for i = 1:length(objects)
    out_points_w = objects(i).out_points_w;
    objects(i).type = 9;
    objects(i).points = out_points_w ./ repmat(sqrt(sum(out_points_w.^2,2)),1,3);
    objects(i).x_w(1:3) = out_points_w(3,:);
    objects(i).x_w(4) = norm(out_points_w(2,:)-out_points_w(3,:));
    objects(i).x_w(5) = norm(out_points_w(6,:)-out_points_w(2,:));
    objects(i).x_w(6) = norm(out_points_w(4,:)-out_points_w(3,:));
    objects(i).x_w(7) = main_angle(views(1)); 
end

end

