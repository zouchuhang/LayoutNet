function [ panoList ] = edgeFromImg2Pano( edges )
%EDGEFROMIMG2PANO Summary of this function goes here
%   Detailed explanation goes here
edgeList = edges.edgeLst;
if size(edgeList,1) == 0
    panoList = [];
    return;
end

vx = edges.vx;
vy = edges.vy;
fov = edges.fov;
[imH, imW, ~] = size(edges.img);

R = (imW/2) / tan(fov/2);

% im is the tangent plane, contacting with ball at [x0 y0 z0]
x0 = R * cos(vy) * sin(vx);
y0 = R * cos(vy) * cos(vx);
z0 = R * sin(vy);
vecposX = [cos(vx) -sin(vx) 0];
vecposY = cross([x0 y0 z0], vecposX);
vecposY = vecposY ./ sqrt(vecposY*vecposY');
Xc = (0 + imW-1)/2;
Yc = (0 + imH-1)/2;

vecx1 = edgeList(:,1)-Xc;
vecy1 = edgeList(:,2)-Yc;
vecx2 = edgeList(:,3)-Xc;
vecy2 = edgeList(:,4)-Yc;

vec1 = repmat(vecx1, [1 3]).*repmat(vecposX, [size(vecx1,1) 1]) ...
     + repmat(vecy1, [1 3]).*repmat(vecposY, [size(vecy1,1) 1]);
vec2 = repmat(vecx2, [1 3]).*repmat(vecposX, [size(vecx2,1) 1]) ...
     + repmat(vecy2, [1 3]).*repmat(vecposY, [size(vecy2,1) 1]);
coord1 = repmat([x0 y0 z0], [size(vec1,1) 1]) + vec1;
coord2 = repmat([x0 y0 z0], [size(vec2,1) 1]) + vec2;

normal = cross( coord1, coord2, 2);
n = sqrt( normal(:,1).^2 + normal(:,2).^2 + normal(:,3).^2);
normal = normal ./ repmat(n, [1 3]);

panoList = [normal coord1 coord2 edgeList(:,end)];

%% lines



end

