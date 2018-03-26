function [ oc ] = occlusionTest( seedrect, testrects )
%OCCLUSIONTEST Summary of this function goes here
%   Detailed explanation goes here
max_min_x = max(testrects(:,1), seedrect(1));
max_min_y = max(testrects(:,2), seedrect(2));
max_min_z = max(testrects(:,3), seedrect(3));
min_max_x = min(testrects(:,4), seedrect(4));
min_max_y = min(testrects(:,5), seedrect(5));
min_max_z = min(testrects(:,6), seedrect(6));

oc = zeros(size(testrects,1),1);
valid = max_min_x<=min_max_x & max_min_y<=min_max_y & max_min_z<=min_max_z;
inter = (min_max_x(valid,1)-max_min_x(valid,1)+0.00001) ...
     .* (min_max_y(valid,1)-max_min_y(valid,1)+0.00001) ...
     .* (min_max_z(valid,1)-max_min_z(valid,1)+0.00001);
seedarea = repmat(prod(seedrect(4:6)-seedrect(1:3)+0.00001), sum(valid),1);
testarea = prod(testrects(valid,4:6)-testrects(valid,1:3)+0.00001, 2);
% oc(valid) = inter./(seedarea+testarea-inter);
oc(valid) = inter./min(seedarea, testarea);
% try
%     oc(valid) = inter./min(seedarea, testarea);
% catch
%     fprintf('error\n');
% end

end

