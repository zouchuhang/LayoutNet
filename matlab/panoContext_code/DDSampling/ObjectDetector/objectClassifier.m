function [ scores, valid ] = objectClassifier(room_obj, objects, model, ClassNames )
%OBJECT_CLASSIFIER Summary of this function goes here
%   Detailed explanation goes here
% ClassNames = model.ClassNames;
% assert(length(ClassNames)==28, 'Classifier involve class number less than 28.');
% fprintf('object classifier: %d categories\n', length(ClassNames));

room_points = room_obj.out_points_w;
room_xyz = [min(room_points, [], 1) max(room_points, [],  1)];
total_fea_dim = comp_object_3d_feature();
num_object = length(objects);
features = zeros(num_object, total_fea_dim);
valid = true(num_object,1);

for oid = 1:num_object
    obj_points = objects(oid).out_points_w;
    obj_xyz = [min(obj_points, [], 1) max(obj_points, [], 1)];
    test_min = (room_xyz(1:3) - obj_xyz(1:3));
    norm_min = (room_xyz(1:3) - obj_xyz(1:3))./(obj_xyz(4:6)-obj_xyz(1:3)+0.00001);
    test_max = (obj_xyz(4:6) - room_xyz(4:6));
    norm_max = (obj_xyz(4:6) - room_xyz(4:6))./(obj_xyz(4:6)-obj_xyz(1:3)+0.00001);
    if any(norm_min>0.2) || any(norm_max>0.2) || any(test_min>20) || any(test_max>20)
        valid(oid) = false;
        continue;
    end

    obj_x = objects(oid).x_w;
    obj_ori_size = obj_x(4:6)';
    obj_angle = obj_x(7);
    features(oid,:) = comp_object_3d_feature(obj_xyz, room_xyz, obj_ori_size, obj_angle);
end

[~,scores,~] = predict(model, features);
if ~exist('ClassNames','var') || isempty(ClassNames)
    ClassNames = model.ClassNames;
end
I = zeros(1,length(ClassNames));
for i = 1:length(ClassNames)
    I(i) = str2double(ClassNames{i});
end
scores(:,I) = scores;

end

