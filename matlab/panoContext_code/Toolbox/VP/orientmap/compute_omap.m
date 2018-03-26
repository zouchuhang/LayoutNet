function [omap, OMAP_FACTOR] = compute_omap(lines, vp, imgsize)

% global OMAP_FACTOR; 
% OMAP_FACTOR = 250 / norm(imgsize(1:2));
OMAP_FACTOR = 1;

[omap omapstrict1 omapstrict2] = compute_orientationmap(lines, vp, imgsize(1:2), OMAP_FACTOR);
% omapstrict1, omapstrict2 are more conservative estimates.
% They should be more reliable but have more uncertain regions.
