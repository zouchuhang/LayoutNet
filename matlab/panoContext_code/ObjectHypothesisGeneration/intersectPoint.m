function [ intersect ] = intersectPoint( line1, line2 )
%INTERSECTPOINT Summary of this function goes here
%   line is represented by two end point

A = (line1(:,1)-line1(:,3))./(line1(:,2)-line1(:,4)+0.00001);
B = (line2(:,1)-line2(:,3))./(line2(:,2)-line2(:,4)+0.00001);

Y = (-(line1(:,1)-line2(:,1)) + (A.*line1(:,2)-B.*line2(:,2)))./(A-B+0.00001);
X = A.*Y - A.*line1(:,2) + line1(:,1);

intersect = [X Y];
end

