function [ feature ] = comp_object_3d_feature( obj_xyz, room_xyz, obj_ori_size, obj_angle )
%COMP_OBJECT_3D_FEATURE Summary of this function goes here
%   Detailed explanation goes here
if nargin==0
    feature = 55;
    return;
end

obj_ctr = (obj_xyz(1:3)+obj_xyz(4:6))/2;
obj_size = obj_xyz(4:6) - obj_xyz(1:3);
room_size = room_xyz(4:6) - room_xyz(1:3);

feature = zeros(1,49);

feature(1:3) = obj_size; % object size
feature(4:6) = obj_size./room_size; % relative size

feature(7:9) = obj_ctr - room_xyz(1:3); % center position to min
feature(10:12) = room_xyz(4:6) - obj_ctr; % center position to max
feature(13:15) = feature(7:9)./obj_size; % center position to min wrt object size
feature(16:18) = feature(10:12)./obj_size; % center position to max wrt object size
feature(19:21) = feature(7:9)./room_size; % center position to min wrt room size

feature(22:24) = obj_xyz(1:3)-room_xyz(1:3); % closest position to min
feature(25:27) = room_xyz(4:6)-obj_xyz(4:6); % closest position to max
feature(28:30) = feature(22:24)./obj_size; % closest position to min wrt object size
feature(31:33) = feature(25:27)./obj_size; % closest position to max wrt object size
feature(34:36) = feature(22:24)./room_size; % closest position to min wrt room size

feature(37:39) = obj_size([1 2 3])./obj_size([2 3 1]); % aspect ratio
sort_obj_size = sort(obj_size, 'ascend');
feature(40:42) = sort_obj_size([1 2 3])./sort_obj_size([2 3 1]); % sorted aspect ratio
feature(43:45) = obj_ori_size([1 2 3])./obj_ori_size([2 3 1]); % original aspect ratio
sort_obj_ori_size = sort(obj_ori_size, 'ascend');
feature(46:48) = sort_obj_ori_size([1 2 3])./sort_obj_ori_size([2 3 1]); % sorted original aspect ratio

sort_obj_size = [sort(obj_size(1:2)) obj_size(3)];
feature(49:51) = sort_obj_size([1 2 3])./sort_obj_size([2 3 1]); % xy sorted aspect ratio
sort_obj_ori_size = [sort(obj_ori_size(1:2)) obj_ori_size(3)];
feature(52:54) = sort_obj_ori_size([1 2 3])./sort_obj_ori_size([2 3 1]); % xy sorted aspect ratio

feature(55) = abs(round(obj_angle/(pi/2))*(pi/2)-obj_angle);

end

