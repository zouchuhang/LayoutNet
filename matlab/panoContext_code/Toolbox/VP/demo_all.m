% demo_all.m
% Demo script for the following operations:
%  - extract line segments
%  - estimate vanishing points
%  - compute omap
%  - generate cuboid hypotheses from omap
%  - sample room hypotheses from line segments
%
% Assumptions about the image
% 1. Manhattan, there are three mutually orthogonal vanishing points.
% 2. There is a vanishing point that is close to vertical.
% 3. Principal point is assumed to be at the center of the image, 
%    i.e., cropped images will not work.
%
% Related publications:
% David C. Lee, Abhinav Gupta, Martial Hebert, and Takeo Kanade.
%    "Estimating Spatial Layout of Rooms using Volumetric Reasoning about 
%    Objects and Surfaces." Advances in Neural Information Processing 
%    Systems 24 (NIPS) 2010.
% David C. Lee, Martial Hebert, and Takeo Kanade. "Geometric Reasoning for
%    Single Image Structure Recovery." IEEE Conference on Computer Vision 
%    and Pattern Recognition (CVPR) 2009.
% 
% 05/2011 David C. Lee (dclee@cs.cmu.edu)
%

%%
% img = imread('uiuc261.jpg');
img = imread('test1.jpg');
% img = sepScene(22).img;
%%
addpath('display');
addpath('lineseg');
addpath('vanishingpoint');
addpath('geometry');
addpath('orientmap');
addpath('genobjhyp');
addpath('genroom');

%% Extract line segments
% Note: 'lines' are better for computing vanishing points.
%       'linesmore' are better for computing orientation maps.
if ispc % windows binary (recommended)
    [lines linesmore] = compute_lines(img);
else % peter kovesi line segment detector
    addpath('lineseg/pkline');
    lines = pkline(rgb2gray(img));
end

disp_lines(img, lines);

%% Compute vanishing point and focal length

[vp f] = compute_vp(lines, size(img));

[lines lines_ex] = taglinesvp(vp, lines);
[linesmore linesmore_ex] = taglinesvp(vp, linesmore);

fprintf('vanishing points: [%.1f,%.1f], [%.1f,%.1f], [%.1f,%.1f]\n',...
    vp{1}(1),vp{1}(2), vp{2}(1),vp{2}(2), vp{3}(1),vp{3}(2));
fprintf('focal length: %.1f\n', f);
disp_vanish(img, lines, vp);
disp_vanish(img, lines, vp); axis auto
disp_vanish(img, linesmore, vp);

%% Compute orientation map
% Note: 'lines' are better for computing vanishing points.
%       'linesmore' are better for computing orientation maps.
% [omap, OMAP_FACTOR] = compute_omap(lines, vp, size(img));
% disp_omap(omap, img, 0.6);
[omapmore, OMAP_FACTOR] = compute_omap(linesmore, vp, size(img));
disp_omap(omapmore, img, 0.6);

%% Generate cuboid hypotheses from omap
cuboidhyp_omap = generate_cuboid_from_omap(omapmore, vp, OMAP_FACTOR);

disp_cubes(cuboidhyp_omap, img, 1); % display all
% disp_cubes(cuboidhyp_omap, img, 2); % quickly flash through
% disp_cubes(cuboidhyp_omap, img, 3); % examine one by one

%% Generate room hypotheses by sampling line segments
roomhyp = sample_roomhyp(1000, linesmore_ex, vp, size(img));

% disp_room(roomhyp, img, 1); % display all
% disp_room(roomhyp, img, 2); % quickly flash through
% disp_room(roomhyp, img, 3); % examine one by one
disp_room(roomhyp(randsample(length(roomhyp),10)), img, 1); % display some

%% Some extra functions that may be useful :)
RUN_EXTRA = 0;
if RUN_EXTRA
    addpath('evalhyp');
    
    cid = 1;
    [omapcube rmapcube] = cube_to_orientmap(cuboidhyp_omap(cid), size(img), vp, OMAP_FACTOR);
    disp_omap(omapcube, img, 0.6);
%     hold on; disp_cubes(cuboidhyp_omap(cid), [], 1);

    rid = 1;
    [omapbld rmapbld] = get_labelimg_frombox(roomhyp(rid).box, size(img), vp, OMAP_FACTOR);
    disp_omap(omapbld, img, 0.6);
%     hold on; disp_room(roomhyp(rid), [], 1);
end

