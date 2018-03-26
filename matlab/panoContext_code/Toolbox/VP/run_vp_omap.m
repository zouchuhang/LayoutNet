% run_vp_omap.m
% Demo script for computing vanishing points and orientation map
% 05/2011 David C. Lee (dclee@cs.cmu.edu)
%
% Assumptions about the image
% 1. Manhattan, there are three mutually orthogonal vanishing points.
% 2. There is a vanishing point that is close to vertical.
% 3. Principal point is assumed to be at the center of the image, 
%    i.e., cropped images will not work.
%

%%
% img = imread('uiuc261.jpg');
img = sepScene(1).img;
%%
addpath('display');
addpath('lineseg');
addpath('vanishingpoint');
addpath('geometry');
addpath('orientmap');

%%
% Note: 'lines' are better for computing vanishing points.
%       'linesmore' are better for computing orientation maps.

if ispc % windows binary (recommended)
    [lines linesmore] = compute_lines(img);
else % peter kovesi line segment detector
    addpath('lineseg/pkline');
    lines = pkline(rgb2gray(img));
end

disp_lines(img, lines);

%%
[vp f] = compute_vp(lines, size(img));

fprintf('vanishing points: [%.1f,%.1f], [%.1f,%.1f], [%.1f,%.1f]\n',...
    vp{1}(1),vp{1}(2), vp{2}(1),vp{2}(2), vp{3}(1),vp{3}(2));
fprintf('focal length: %.1f\n', f);

lines_orig = lines; 
[lines lines_ex] = taglinesvp(vp, lines_orig);
linesmore_orig = linesmore; 
[linesmore linesmore_ex] = taglinesvp(vp, linesmore_orig);

%%
disp_vanish(img, lines, vp);
disp_vanish(img, lines, vp); axis auto
disp_vanish(img, linesmore, vp);

%%
% [omap, OMAP_FACTOR] = compute_omap(lines, vp, size(img));
% disp_omap(omap, img, 0.6);

[omapmore, OMAP_FACTOR] = compute_omap(linesmore, vp, size(img));
disp_omap(omapmore, img, 0.6);
