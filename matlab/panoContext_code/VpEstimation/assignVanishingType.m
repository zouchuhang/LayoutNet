function [ type, typeCost ] = assignVanishingType( lines, vp, tol, area )
%ASSIGNVANISHINGTYPE Summary of this function goes here
%   Detailed explanation goes here
if nargin<=3
    area = 10;
end

numLine = size(lines,1);
numVP = size(vp, 1);
typeCost = zeros(numLine, numVP);
% perpendicular 
for vid = 1:numVP
    cosint = dot( lines(:,1:3), repmat( vp(vid,:), [numLine 1]), 2);
    typeCost(:,vid) = asin(abs(cosint));
end
% infinity
for vid = 1:numVP
    valid = true(numLine,1);
    for i = 1:numLine
        us = lines(i,5);
        ue = lines(i,6);
        u = [us;ue]*2*pi-pi;
        v = computeUVN(lines(i,1:3), u, lines(i,4));
        xyz = uv2xyzN([u v], lines(i,4));
        x = linspace(xyz(1,1),xyz(2,1),100);
        y = linspace(xyz(1,2),xyz(2,2),100);
        z = linspace(xyz(1,3),xyz(2,3),100);
        xyz = [x' y' z'];
        xyz = xyz ./ repmat(sqrt(sum(xyz.^2,2)),[1 3]);
        ang = acos( abs(dot(xyz, repmat(vp(vid,:), [100 1]), 2)));
        valid(i) = ~any(ang<area*pi/180);
    end
    typeCost(~valid,vid) = 100;
end

[I, type] = min(typeCost,[],2);
type(I>tol) = numVP+1;

end

