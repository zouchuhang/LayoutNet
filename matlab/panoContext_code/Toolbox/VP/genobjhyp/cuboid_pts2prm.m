function prm = cuboid_pts2prm(pts, vp)
% converts eight 2D points to five parameters of cuboid
% (performs nonlinear optimization)
%
% prm: [x, y, w, h, d]
% x,y: 2D coordinates of corner1 (top,left,front corner)
% w: length (in pixels) from corner1 to corner2
% h: length (in pixels) from corner1 to corner3
% d: length (in pixels) from corner1 to corner5
%
% pts: [8 x 2] 2D coordinates of eight points on cuboid
%

prm_init = zeros(1,5);
prm_init(1:2) = pts(1,:);
prm_init(3) = norm(pts(1,:) - pts(2,:));
prm_init(4) = norm(pts(1,:) - pts(3,:));
prm_init(5) = norm(pts(1,:) - pts(5,:));

%%
prm = lsqnonlin(@(x) ptserr(x, pts, vp), prm_init);


%%
function err = ptserr(prm, pts, vp)
% mean squared error of pts & cuboid encoded by prm

pts1 = cuboid_prm2pts(prm, vp);
err = sum((pts1 - pts).^2, 2);

