function prm = cuboid_pts2prm_incomp(cuboidhyp, vp)
% given incomplete number of points defining a cuboid
% converts to five parameters of cuboid
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

%% initial parameters
pts = cat(1, cuboidhyp.junc3.pt);

prm_init = zeros(1,5);
prm_init(1:2) = mean(pts, 1); % rough initialization
r = mean(std(pts, 0, 1)); % rough initialization
prm_init(3) = r;
prm_init(4) = r;
prm_init(5) = r;

%%
ptexist = false(1,8); % list of corners that exist
for i = 1:8
    if ~isempty(cuboidhyp.junc3(i).pt)
        ptexist(i) = 1;
    end
end

%%
prm = lsqnonlin(@(x) ptserr(x, pts, vp, ptexist), prm_init);
% prm = lsqnonlin(@(x) ptserr(x, pts, vp, ptexist), prm_init, ...
%     [-Inf -Inf 5 5 5], [Inf Inf Inf Inf Inf]);


%%
function err = ptserr(prm, pts, vp, ptexist)
% mean squared error of pts & cuboid encoded by prm

pts1 = cuboid_prm2pts(prm, vp);
% err = sum((pts1(ptexist,:) - pts).^2, 2);
diff = pts1(ptexist,:) - pts;
err = diff(:);

