function [ objects ] = floorRectangle2hypothesis( xyzBox, rule, room_points )
%RECTANGLE2HYPOTHESIS Summary of this function goes here
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
LRDU_SURF(5,:) = [1 2 3 4];
LRDU_SURF(6,:) = [1 2 4 3];

normals = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
normal_ids = [1 1 2 2 3 3];
main_angle = [pi/2 -pi/2 -pi 0 0];

objects = repmat(struct('out_points_w',[],'name',[],'type',[],'points',[],'x_w',[]),0,1);
object = struct('out_points_w',[],'name','other','type',[],'points',[],'x_w',zeros(7,1));
xyzBox = reshape(xyzBox', 3, 4)';

%% if left two points in lhs wall, add a hypothesis
[inside,~,~] = insideCone( wall_points(:,:,LRDU_SURF(rule,1)), xyzBox([1 4],:), 0 );
if all(inside)   
    LID = LRDU_SURF(rule,1);
    wall_point = wall_points(:,:,LID);
    point5 = LineFaceIntersection(wall_point(1,:), normals(LID,:), [0 0 0], xyzBox(1,:));
    point4 = LineFaceIntersection(wall_point(1,:), normals(LID,:), [0 0 0], xyzBox(4,:));
    point6 = LineFaceIntersection(point5, normals(rule,:), [0 0 0], xyzBox(2,:));
    point1 = LineFaceIntersection(point4, normals(rule,:), [0 0 0], xyzBox(3,:));
    point3 = point4;    point2 = point1;
    point3(normal_ids(rule)) = wall_points(1,normal_ids(rule),rule);
    point2(normal_ids(rule)) = wall_points(1,normal_ids(rule),rule);
    
    object.out_points_w = [point1; point2; point3; point4; point5; point6];
    object.type = 13;
    object.points = object.out_points_w ./ repmat(sqrt(sum(object.out_points_w.^2,2)),1,3);
    object.x_w(1:3) = point3;
    object.x_w(4) = norm(point2-point3);
    object.x_w(5) = norm(point5-point4);
    object.x_w(6) = norm(point4-point3);
    object.x_w(7) = main_angle(rule); 
    objects(end+1,1) = object;
end

%% if right two points in rhs wall, add a hypothesis
[inside,~,~] = insideCone( wall_points(:,:,LRDU_SURF(rule,2)), xyzBox([2 3],:), 0 );
if all(inside)
    RID = LRDU_SURF(rule,2);
    wall_point = wall_points(:,:,RID);
    point6 = LineFaceIntersection(wall_point(1,:), normals(RID,:), [0 0 0], xyzBox(2,:));
    point1 = LineFaceIntersection(wall_point(1,:), normals(RID,:), [0 0 0], xyzBox(3,:));
    point5 = LineFaceIntersection(point6, normals(rule,:), [0 0 0], xyzBox(1,:));
    point4 = LineFaceIntersection(point1, normals(rule,:), [0 0 0], xyzBox(4,:));
    point3 = point4;    point2 = point1;
    point3(normal_ids(rule)) = wall_points(1,normal_ids(rule),rule);
    point2(normal_ids(rule)) = wall_points(1,normal_ids(rule),rule);
    
    object.out_points_w = [point1; point2; point3; point4; point5; point6];
    object.type = 13;
    object.points = object.out_points_w ./ repmat(sqrt(sum(object.out_points_w.^2,2)),1,3);
    object.x_w(1:3) = point3;
    object.x_w(4) = norm(point2-point3);
    object.x_w(5) = norm(point5-point4);
    object.x_w(6) = norm(point4-point3);
    object.x_w(7) = main_angle(rule); 
    objects(end+1,1) = object;
end

%% if bottom two points in bottom wall, add a hypothesis
[inside,~,~] = insideCone( wall_points(:,:,LRDU_SURF(rule,3)), xyzBox([4 3],:), 0 );
if all(inside)
    BID = LRDU_SURF(rule,3);
    wall_point = wall_points(:,:,BID);
    point4 = LineFaceIntersection(wall_point(1,:), normals(BID,:), [0 0 0], xyzBox(4,:));
    point1 = LineFaceIntersection(wall_point(1,:), normals(BID,:), [0 0 0], xyzBox(3,:));
    point5 = LineFaceIntersection(point4, normals(rule,:), [0 0 0], xyzBox(1,:));
    point6 = LineFaceIntersection(point1, normals(rule,:), [0 0 0], xyzBox(2,:));
    point3 = point4;    point2 = point1;
    point3(normal_ids(rule)) = wall_points(1,normal_ids(rule),rule);
    point2(normal_ids(rule)) = wall_points(1,normal_ids(rule),rule);
    
    object.out_points_w = [point1; point2; point3; point4; point5; point6];
    object.type = 13;
    object.points = object.out_points_w ./ repmat(sqrt(sum(object.out_points_w.^2,2)),1,3);
    object.x_w(1:3) = point3;
    object.x_w(4) = norm(point2-point3);
    object.x_w(5) = norm(point5-point4);
    object.x_w(6) = norm(point4-point3);
    object.x_w(7) = main_angle(rule); 
    objects(end+1,1) = object;
end

%% if top two points in upper wall, add a hypothesis
[inside,~,~] = insideCone( wall_points(:,:,LRDU_SURF(rule,4)), xyzBox([1 2],:), 0 );
if all(inside)
    UID = LRDU_SURF(rule,4);
    wall_point = wall_points(:,:,UID);
    point5 = LineFaceIntersection(wall_point(1,:), normals(UID,:), [0 0 0], xyzBox(1,:));
    point6 = LineFaceIntersection(wall_point(1,:), normals(UID,:), [0 0 0], xyzBox(2,:));
    point4 = LineFaceIntersection(point5, normals(rule,:), [0 0 0], xyzBox(4,:));
    point1 = LineFaceIntersection(point6, normals(rule,:), [0 0 0], xyzBox(3,:));
    point3 = point4;    point2 = point1;
    point3(normal_ids(rule)) = wall_points(1,normal_ids(rule),rule);
    point2(normal_ids(rule)) = wall_points(1,normal_ids(rule),rule);
    
    object.out_points_w = [point1; point2; point3; point4; point5; point6];
    object.type = 13;
    object.points = object.out_points_w ./ repmat(sqrt(sum(object.out_points_w.^2,2)),1,3);
    object.x_w(1:3) = point3;
    object.x_w(4) = norm(point2-point3);
    object.x_w(5) = norm(point5-point4);
    object.x_w(6) = norm(point4-point3);
    object.x_w(7) = main_angle(rule); 
    objects(end+1,1) = object;
end

end

