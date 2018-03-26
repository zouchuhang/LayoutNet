function pts = cuboid_prm2pts(prm, vp)
% converts five parameters of cuboid to eight 2D points
%
% prm: [x, y, w, h, d]
% x,y: 2D coordinates of corner1 (top,left,front corner)
% w: length (in pixels) from corner1 to corner2
% h: length (in pixels) from corner1 to corner3
% d: length (in pixels) from corner1 to corner5
%
% pts: [8 x 2] 2D coordinates of eight points on cuboid
%   5---6   1------2
%  /   /|   |      |
% 1---2 8   3------4    ....
% |   |/     \    /
% 3---4       7--8
%
% 1: front,top,left
% 2: front,top,right
% ...
% 8: back,bottom,right

pts = zeros(8,2);

pts(1,:) = [prm(1) prm(2)];
pts(2,:) = towardsright(pts(1,:), vp, prm(3));
pts(3,:) = towardsdown(pts(1,:), vp, prm(4));
pts(5,:) = towardsback(pts(1,:), vp, prm(5));
pts(4,:) = line_intersect(pts(2,:), vp{1}, pts(3,:), vp{2});
pts(6,:) = line_intersect(pts(2,:), vp{3}, pts(5,:), vp{2});
pts(7,:) = line_intersect(pts(3,:), vp{3}, pts(5,:), vp{1});
pts(8,:) = mean(...
    [line_intersect(pts(4,:), vp{3}, pts(6,:), vp{1}); ...
     line_intersect(pts(7,:), vp{2}, pts(6,:), vp{1}); ...
     line_intersect(pts(4,:), vp{3}, pts(7,:), vp{2})], 1);

