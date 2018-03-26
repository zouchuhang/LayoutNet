function [vp f] = vanish_from_minevidence(lines, imgsize)
%
% find three orthogonal vanishing points by 
% 1. findinhg vertical vp.
% 2. selecting one additional line.
% 3. pick the best f for the current triplet of vanishing points.
% 4. go to 2 and pick the best additional line.

%%
THRES_THETA = 10;
THRES_THETA_GUESS = 20;

%%
for i = 1:length(lines)
    lines(i).length = norm(lines(i).point1 - lines(i).point2);
end
global maxlinelength;
maxlinelength = max([lines.length]);

%% find vertical vanishing point
guessvp = [0 -1e5]; % initial guess of the vertical vp
lc1_guess = line_belongto_vp(lines, guessvp, THRES_THETA_GUESS);

vp1 = find_single_vp(lines(lc1_guess), THRES_THETA);
lc1 = line_belongto_vp(lines, vp1, THRES_THETA);

% global img; disp_lines(img, lines);
% global img; disp_lines(img, lines(lc1_guess));
% global img; disp_lines(img, lines(lc1)); plot(vp1(1), vp1(2), 'x', 'LineWidth', 5, 'MarkerSize', 20);

%% find the other two vanishing points by 
[vp2 vp3 f] = find_other_vp(lines(~lc1), vp1, THRES_THETA, imgsize);

%%
INFDIST = 1e7; % distance to a far away point...
if ~isreal(vp1), vp1 = INFDIST * vp1 / j; end
if ~isreal(vp2), vp2 = INFDIST * vp2 / j; end
if ~isreal(vp3), vp3 = INFDIST * vp3 / j; end

%% order vp
imgcenter = [imgsize(2) imgsize(1)]/2;
if norm(vp2-imgcenter) < norm(vp3-imgcenter)
    temp = vp2;
    vp2 = vp3;
    vp3 = temp;
end

vp = {vp1; vp2; vp3};


%%
function [vp2 vp3 f] = find_other_vp(lines_other, vp1, THRES_THETA, imgsize)

NUMITER = 200; % number of lines to sample
FSWEEP = 50:50:1000; % focal length to test

if length(lines_other) > NUMITER
    linesample = randsample(length(lines_other), NUMITER);
else
    linesample = 1:length(lines_other);
end

scorearray = zeros(length(linesample), length(FSWEEP));
for i = 1:length(linesample)
for j = 1:length(FSWEEP)
    line = lines_other(linesample(i));
    f = FSWEEP(j);
    
    [vp2 vp3] = vpset_from_vp1_line_f(vp1, line, f, imgsize);

    % line assignment...
    [vv2 vv2angle] = line_belongto_vp(lines_other, vp2, THRES_THETA);
    [vv3 vv3angle] = line_belongto_vp(lines_other, vp3, THRES_THETA);
    vvboth = vv2 & vv3;
    vvboth2 = vvboth & (vv2angle<vv3angle);
    vvboth3 = vvboth & (vv2angle>=vv3angle);
    vvvv2 = (vv2 & ~vv3) | vvboth2;
    vvvv3 = (vv3 & ~vv2) | vvboth3;

    % score
    score = getvote(vp2, lines_other(vvvv2), THRES_THETA) + ...
            getvote(vp3, lines_other(vvvv3), THRES_THETA);
    
    scorearray(i,j) = score;
end
end

[c idx] = max(scorearray(:));
[row col] = ind2sub(size(scorearray), idx);

line = lines_other(linesample(row));
f = FSWEEP(col);
[vp2 vp3] = vpset_from_vp1_line_f(vp1, line, f, imgsize);


%%
function [vp2 vp3] = vpset_from_vp1_line_f(vp1, line, f, imgsize)

[vl p1 p2] = vp2vl(vp1, f, imgsize);

% intersection of chosen line and vl
% [intpt degen] = line_intersect(p1, p2, line.point1, line.point2);
intpt = line_intersect_forvp(p1, p2, line.point1, line.point2);

% second vp is the intpt
% third vp is orthogonal to the 1st and 2nd
vp2 = intpt;
% Is3 = cross([vp1(:)-imgsize([2 1])'/2; f], [vp2(:)-imgsize([2 1])'/2; f]);
Is3 = cross( xy2is(vp1, f, imgsize), xy2is(vp2, f, imgsize) );
vp3 = Is3(1:2)' / Is3(3) * f + imgsize([2 1])/2;
    

%%
function [vl p1 p2] = vp2vl(vp, f, imgsize)
% returns the equation for the vanishing line of planes
% that are orthogonal to the lines which are
% associated with the vanishing point.
%
% vl [3x1]: vl(1)*x + vl(2)*y + vl(3) = 0
% p1, p2: two points on line

% principal point (assume at image center)
cx = imgsize(2)/2;
cy = imgsize(1)/2;

% if vp is at infinity...
if ~isreal(vp)
%     vl is a line passing through image center.
    a = vp(1) / j;
    b = vp(2) / j;
    c = - a*cx - b*cy;
    vl = [a b c];
    p1 = [cx cy];
    p2 = [cx+b cy-a];
    return;
end

% normalized coord of vp
nx = vp(1) - cx;
ny = vp(2) - cy;

% ray
% r = [nx ny f]';

% ortho plane
% p = r / norm(r);

% direction of intersecting line
% [nx ny f]' x [0 0 1]' = [ny -nx 0]

% one point on line
% on imaging plane: z0 = f
% on ortho plane:   nx*nx0 + ny*ny0 + f*z0 = 0
if abs(nx) > abs(ny)
    z0 = f;
    ny0 = 0;
    nx0 = - (ny * ny0 + f * f) / nx;
else
    z0 = f;
    nx0 = 0;
    ny0 = - (nx * nx0 + f * f) / ny;
end
x0 = nx0 + cx;
y0 = ny0 + cy;

% line eq in 3D
% [nx0 ny0 z0]' + lambda * [ny -nx 0]

% line eq in 2D
% [x0 y0]' + lambda * [ny -nx]
% <==>
% ax + by + c = 0;
a = nx;
b = ny;
c = - a * x0 - b * y0;

vl = [a b c]';

p1 = [x0 y0];
p2 = [x0 y0] + [ny -nx];


%%
function is = xy2is(xy, f, imgsize)

if ~isreal(xy)
    is = [xy(1)/j; xy(2)/j; 0];
    is = is / norm(is);
else
    is = [xy(1)-imgsize(2)/2; xy(2)-imgsize(1)/2; f];
    is = is / norm(is);
end






