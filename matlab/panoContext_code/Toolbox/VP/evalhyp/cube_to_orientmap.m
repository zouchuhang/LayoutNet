function [omap rmap] = cube_to_orientmap(cube, imgsize, vp, OMAP_FACTOR)


% % rmap id (cube surfaceid)
% botid = 1; % 1: bottom
% frtid = 2; % 2: front
% lftid = 3; % 3: left
% rhtid = 4; % 4: right
% bckid = 5; % 5: back
% topid = 6; % 6: top


%% resize cube, imgsize, vp
for i = 1:8
    cube.junc3(i).pt = cube.junc3(i).pt * OMAP_FACTOR;
end
vp{1} = vp{1} * OMAP_FACTOR;
vp{2} = vp{2} * OMAP_FACTOR;
vp{3} = vp{3} * OMAP_FACTOR;
imgsize = ceil(imgsize * OMAP_FACTOR);

%%

omap = zeros(imgsize(1), imgsize(2), 3);
rmap = zeros(imgsize(1), imgsize(2), 1);

%% determine visible surfaces
vissurf = get_cube_visible_surface(cube, vp);
% cube surfaceid
% botid = 1; % 1: bottom
% frtid = 2; % 2: front
% lftid = 3; % 3: left
% rhtid = 4; % 4: right
% bckid = 5; % 5: back
% topid = 6; % 6: top


%% gen map for visible surfaces
% [1234(front), 1256(top), 3478(bot), 1357(left), 2468(right)]
% front 1,2,3,4 -- always visible
if vissurf(2) == 1
    pt = cat(1, cube.junc3([1 2 4 3]).pt);
    mask = poly2mask(pt(:,1), pt(:,2), imgsize(1), imgsize(2));
    omap(:,:,3) = omap(:,:,3) | mask;
    rmap(mask) = 2;
end
% top 1,2,5,6
if vissurf(6) == 1
    pt = cat(1, cube.junc3([1 2 6 5]).pt);
    mask = poly2mask(pt(:,1), pt(:,2), imgsize(1), imgsize(2));
    omap(:,:,1) = omap(:,:,1) | mask;
    rmap(mask) = 6;
end
% bot 3,4,7,8
if vissurf(1) == 1
    pt = cat(1, cube.junc3([3 4 8 7]).pt);
    mask = poly2mask(pt(:,1), pt(:,2), imgsize(1), imgsize(2));
    omap(:,:,1) = omap(:,:,1) | mask;
    rmap(mask) = 1;
end
% left 1,3,5,7
if vissurf(3) == 1
    pt = cat(1, cube.junc3([1 3 7 5]).pt);
    mask = poly2mask(pt(:,1), pt(:,2), imgsize(1), imgsize(2));
    omap(:,:,2) = omap(:,:,2) | mask;
    rmap(mask) = 3;
end
% right 2,4,6,8
if vissurf(4) == 1
    pt = cat(1, cube.junc3([2 4 8 6]).pt);
    mask = poly2mask(pt(:,1), pt(:,2), imgsize(1), imgsize(2));
    omap(:,:,2) = omap(:,:,2) | mask;
    rmap(mask) = 4;
end


