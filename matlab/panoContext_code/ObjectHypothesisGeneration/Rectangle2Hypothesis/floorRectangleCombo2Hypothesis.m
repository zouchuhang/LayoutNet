function [ objects ] = floorRectangleCombo2Hypothesis( xyzBoxes, views, room_points )
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

if views(2)==5
    t = views(2);
    views(2) = views(1);
    views(1) = t;
    t = xyzBoxes(1:4,:);
    xyzBoxes(1:4,:) = xyzBoxes(5:8,:);
    xyzBoxes(5:8,:) = t;
end

xyzBoxes_norm = xyzBoxes./repmat(sqrt(sum(xyzBoxes.^2,2)),1,3);
uvBoxes = xyz2uvN(xyzBoxes_norm);

rules(1).pid = [4 1 3 2 8 7; 4 1 5 6 8 7];
rules(2).pid = [2 3 1 4 8 7; 2 3 5 6 8 7];
rules(3).pid = [3 4 2 1 8 7; 3 4 5 6 8 7];
rules(4).pid = [1 2 4 3 8 7; 1 2 5 6 8 7];

pid = rules(views(2)).pid;
%% on floor?
[inside,~,~] = insideCone( wall_points(:,:,5), xyzBoxes, 0 );
if any(inside)
    point3 = LineFaceIntersection(uni_room_points(1,:), [0 0 1], [0 0 0], xyzBoxes(pid(2,5),:));
    point2 = LineFaceIntersection(uni_room_points(1,:), [0 0 1], [0 0 0], xyzBoxes(pid(2,6),:));
    point4_1 = LineFaceIntersection(point3, normals(views(2),:), [0 0 0], xyzBoxes(pid(1,3),:));
    point1_1 = LineFaceIntersection(point2, normals(views(2),:), [0 0 0], xyzBoxes(pid(1,4),:));
    point4_2 = LineFaceIntersection(point3, normals(views(2),:), [0 0 0], xyzBoxes(pid(2,3),:));
    point1_2 = LineFaceIntersection(point2, normals(views(2),:), [0 0 0], xyzBoxes(pid(2,4),:));
    
    point4 = point3;    point1 = point2;
    val = (point4_1(3) + point1_1(3) + point4_2(3) + point1_2(3))/4;
    point4(3) = val;    point1(3) = val;
    
    point5 = LineFaceIntersection(point4, [0 0 1], [0 0 0], xyzBoxes(pid(2,1),:));
    point6 = LineFaceIntersection(point1, [0 0 1], [0 0 0], xyzBoxes(pid(2,2),:));
    
    id = normal_ids(views(2));
    val = (point5(id) + point6(id))/2;
    point5 = point4;    point6 = point1;
    point5(id) = val;   point6(id) = val;
    
    object.out_points_w = [point1; point2; point3; point4; point5; point6];
    objects(end+1,1) = object;
end

%% on middle wall
[inside,~,~] = insideCone( wall_points(:,:,views(2)), xyzBoxes, 0 );
if any(inside)
    wall_point = wall_points(:,:,views(2));
    point5 = LineFaceIntersection(wall_point(1,:), normals(views(2),:), [0 0 0], xyzBoxes(pid(2,1),:));
    point6 = LineFaceIntersection(wall_point(1,:), normals(views(2),:), [0 0 0], xyzBoxes(pid(2,2),:));
    
    point4_1 = LineFaceIntersection(point5, [0 0 1], [0 0 0], xyzBoxes(pid(1,3),:));
    point1_1 = LineFaceIntersection(point6, [0 0 1], [0 0 0], xyzBoxes(pid(1,4),:));
    point4_2 = LineFaceIntersection(point5, [0 0 1], [0 0 0], xyzBoxes(pid(2,3),:));
    point1_2 = LineFaceIntersection(point6, [0 0 1], [0 0 0], xyzBoxes(pid(2,4),:));
    
    id = normal_ids(views(2));
    val = (point4_1(id) + point1_1(id) + point4_2(id) + point1_2(id))/4;
    point4 = point5;    point1 = point6;
    point4(id) = val;   point1(id) = val;
    
    point3 = LineFaceIntersection(point4, normals(views(2),:), [0 0 0], xyzBoxes(pid(2,5),:));
    point2 = LineFaceIntersection(point1, normals(views(2),:), [0 0 0], xyzBoxes(pid(2,6),:));
    val = (point3(3) + point2(3))/2;
    point3 = point4;    point2 = point1;
    point3(3) = val;    point2(3) = val;
    
    object.out_points_w = [point1; point2; point3; point4; point5; point6];
    objects(end+1,1) = object;
end

%% on lhs wall of the lhs views
[inside,~,~] = insideCone( wall_points(:,:,LRDU_SURF(views(2),1)), xyzBoxes, 0 );
if any(inside)
    LID = LRDU_SURF(views(2),1);
    wall_point = wall_points(:,:,LID);
    
    point5 = LineFaceIntersection(wall_point(1,:), normals(LID,:), [0 0 0], xyzBoxes(pid(2,1),:));
    point4_1 = LineFaceIntersection(wall_point(1,:), normals(LID,:), [0 0 0], xyzBoxes(pid(1,3),:));
    point4_2 = LineFaceIntersection(wall_point(1,:), normals(LID,:), [0 0 0], xyzBoxes(pid(2,3),:));
    point3 = LineFaceIntersection(wall_point(1,:), normals(LID,:), [0 0 0], xyzBoxes(pid(2,5),:));
    
    point4 = (point4_1 + point4_2)/2;
    d3 = 3;
    d5 = setdiff([1 2], normal_ids(LID));
    t = point5(d5); point5 = point4; point5(d5) = t;
    t = point3(d3); point3 = point4; point3(d3) = t;
    
    point6 = LineFaceIntersection(point5, [0 0 1], [0 0 0], xyzBoxes(pid(2,2),:));
    point2 = LineFaceIntersection(point3, normals(views(2),:), [0 0 0], xyzBoxes(pid(2,6),:));
    id = setdiff([1 2], normal_ids(views(2)));
    val = (point6(id) + point2(id))/2;
    point6 = point5;    point6(id) = val;
    point1 = point4;    point1(id) = val;
    point2 = point3;    point2(id) = val;
    
    object.out_points_w = [point1; point2; point3; point4; point5; point6];
    objects(end+1,1) = object; 
     
end

%% on rhs wall of the rhs views
[inside,~,~] = insideCone( wall_points(:,:,LRDU_SURF(views(2),2)), xyzBoxes, 0 );
if any(inside)
    RID = LRDU_SURF(views(2),2);
    wall_point = wall_points(:,:,RID);
    
    point6 = LineFaceIntersection(wall_point(1,:), normals(RID,:), [0 0 0], xyzBoxes(pid(2,2),:));
    point1_1 = LineFaceIntersection(wall_point(1,:), normals(RID,:), [0 0 0], xyzBoxes(pid(1,4),:));
    point1_2 = LineFaceIntersection(wall_point(1,:), normals(RID,:), [0 0 0], xyzBoxes(pid(2,4),:));
    point2 = LineFaceIntersection(wall_point(1,:), normals(RID,:), [0 0 0], xyzBoxes(pid(2,6),:));
    
    point1 = (point1_1 + point1_2)/2;
    d2 = 3;
    d6 = setdiff([1 2], normal_ids(RID));
    t = point6(d6); point6 = point1; point6(d6) = t;
    t = point2(d2); point2 = point1; point2(d2) = t;
    
    point5 = LineFaceIntersection(point6, [0 0 1], [0 0 0], xyzBoxes(pid(2,1),:));
    point3 = LineFaceIntersection(point2, normals(views(2),:), [0 0 0], xyzBoxes(pid(2,5),:));
    id = setdiff([1 2], normal_ids(views(2)));
    val = (point5(id) + point3(id))/2;
    point5 = point6;    point5(id) = val;
    point4 = point1;    point4(id) = val;
    point3 = point2;    point3(id) = val;
    
    object.out_points_w = [point1; point2; point3; point4; point5; point6];
    objects(end+1,1) = object; 
  
end

for i = 1:length(objects)
    out_points_w = objects(i).out_points_w;
    objects(i).type = 13;
    objects(i).points = out_points_w ./ repmat(sqrt(sum(out_points_w.^2,2)),1,3);
    objects(i).x_w(1:3) = out_points_w(3,:);
    objects(i).x_w(4) = norm(out_points_w(2,:)-out_points_w(3,:));
    objects(i).x_w(5) = norm(out_points_w(6,:)-out_points_w(2,:));
    objects(i).x_w(6) = norm(out_points_w(4,:)-out_points_w(3,:));
    objects(i).x_w(7) = main_angle(views(2)); 
end

end

