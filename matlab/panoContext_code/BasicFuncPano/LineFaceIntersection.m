function [ X ] = LineFaceIntersection( faceX, faceN, lineX, lineD )
%LINEFACEINTERSECTION Intersection of a plane and a line
%   faceX: a point on the plane
%   faceN: the normal direction of the plane
%   lineX: a point on the line
%   lineD: direction of the line
%

% A = dot(faceN,lineD,2)/lineD(1);
% B = -dot(faceX,faceN,2) + (dot(lineX(2:3),faceN(2:3),2)) - lineX(1)/lineD(1)*(dot(lineD(2:3),faceN(2:3),2));
A = sum(faceN.*lineD)/lineD(1);
B = -sum(faceX.*faceN) + (sum(lineX(2:3).*faceN(2:3),2)) - lineX(1)/lineD(1)*(sum(lineD(2:3).*faceN(2:3),2));

x = -B/A;

y = (x-lineX(1))/lineD(1)*lineD(2) + lineX(2);
z = (x-lineX(1))/lineD(1)*lineD(3) + lineX(3);
X = [x y z];
end

