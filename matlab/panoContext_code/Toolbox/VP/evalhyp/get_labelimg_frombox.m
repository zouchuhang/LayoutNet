function [orientlabel_img regionlabel_img] = get_labelimg_frombox(box, imgsize, vp, OMAP_FACTOR)
% orientlabel_img:
% 1 for floor and ceiling
% 2 for walls facing left/right
% 3 for walls facing front/back
%
% regionlabel_img:
% 1~length(box) for walls
% length(box)+1 for floor
% length(box)+2 for ceiling

%% resize with OMAP_FACTOR
imgwidth = imgsize(2);
imgheight = imgsize(1);

if exist('OMAP_FACTOR', 'var')
    for i = 1:length(box)
        box(i).p1 = box(i).p1 * OMAP_FACTOR;
        box(i).p2 = box(i).p2 * OMAP_FACTOR;
        box(i).p3 = box(i).p3 * OMAP_FACTOR;
        box(i).p4 = box(i).p4 * OMAP_FACTOR;
    end
    imgwidth = ceil(imgwidth * OMAP_FACTOR);
    imgheight = ceil(imgheight * OMAP_FACTOR);
    vp{1} = vp{1} * OMAP_FACTOR;
    vp{2} = vp{2} * OMAP_FACTOR;
    vp{3} = vp{3} * OMAP_FACTOR;
end

%% initialize orientlabel_img
% 1 for floor and ceiling
% 2 for walls facing left/right
% 3 for walls facing front/back
orientlabel_img = zeros(imgheight, imgwidth);
orientlabel_img = orientlabel_img + 1; % default label for floor/ceiling


%% initialize regionlabel_img
% 1~length(box) for walls
% length(box)+1 for floor
% length(box)+2 for ceiling

if nargout==2
    [p1ext p2ext] = extline(vp{2}, vp{3}, imgwidth, imgheight);
    if p1ext(1) > p2ext(1), aaa = p1ext; p1ext = p2ext; p2ext = aaa; end

    regionlabel_img = poly2mask([p1ext(1) p2ext(1) imgwidth+20 1-20], ...
        [p1ext(2) p2ext(2) -20 -20], ...
        imgheight, imgwidth);

    regionlabel_img = regionlabel_img + length(box) + 1;
end

%%
BUFF = 20;
for j = 1:length(box)
    if j==1 % leftmost box
        poly = cat(1, [1-BUFF 1-BUFF], [1-BUFF imgheight+BUFF], ...
            box(j).p2, box(j).p4, box(j).p3, box(j).p1);
        % ***HACK*** gets rid of pesky black space on border -- dclee 11/22/2010
        poly(poly==1) = 0;
        poly(poly==imgheight) = imgheight+1;
        poly(poly==imgwidth) = imgwidth+1;
        % ***HACK
        BW = poly2mask(poly(:,1), poly(:,2), imgheight, imgwidth);
%         BW = poly2mask([1-20; 1-20; box(j).p2(1); box(j).p4(1); box(j).p3(1); box(j).p1(1)], ...
%             [1-20; imgheight+20; box(j).p2(2); box(j).p4(2); box(j).p3(2); box(j).p1(2)], ...
%             imgheight, imgwidth);
    elseif j==length(box) % rightmost box
        poly = cat(1, [imgwidth+BUFF imgheight+BUFF], [imgwidth+BUFF 1-BUFF], ...
            box(j).p3, box(j).p1, box(j).p2, box(j).p4);
        % ***HACK*** gets rid of pesky black space on border -- dclee 11/22/2010
        poly(poly==1) = 0;
        poly(poly==imgheight) = imgheight+1;
        poly(poly==imgwidth) = imgwidth+1;
        % ***HACK
        BW = poly2mask(poly(:,1), poly(:,2), imgheight, imgwidth);
%         BW = poly2mask([box(j).p1(1); box(j).p2(1); box(j).p4(1); imgwidth+20; imgwidth+20; box(j).p3(1)], ...
%             [box(j).p1(2); box(j).p2(2); box(j).p4(2); imgheight+20; 1-20; box(j).p3(2)], ...
%             imgheight, imgwidth);
    else
        poly = cat(1, ...
            box(j).p3, box(j).p1, box(j).p2, box(j).p4);
        BW = poly2mask(poly(:,1), poly(:,2), imgheight, imgwidth);
%         BW = poly2mask([box(j).p1(1); box(j).p2(1); box(j).p4(1);  box(j).p3(1)], ...
%             [box(j).p1(2); box(j).p2(2); box(j).p4(2);  box(j).p3(2)], ...
%             imgheight, imgwidth);
    end
    if nargout==2
    	regionlabel_img(BW(:)) = j;
    end
	orientlabel_img(BW(:)) = box(j).orient;
end

% orientlabel_img = orientlabel_img(:);
% regionlabel_img = regionlabel_img(:);

orientlabel_img = cat(3, orientlabel_img==1, orientlabel_img==2, orientlabel_img==3);

%%
function [p1ext p2ext] = extline(p1, p2, imgwidth, imgheight)
% p1: [x y]
% p2: [x y]

dir = p2 - p1;

% intersection with top
% p1(1) + alpha*dir(2) = 1
% int1 = p1 + alpha*dir
alpha = (1 - p1(2)) / dir(2);
intp{1} = p1 + alpha*dir;
intp{1}(2) = 1; % numerical precision issues....

% intersection with bottom
% p1(1) + alpha*dir(2) = imgheight
% int2 = p1 + alpha*dir
alpha = (imgheight - p1(2)) / dir(2);
intp{2} = p1 + alpha*dir;
intp{2}(2) = imgheight; % numerical precision issues....

% intersection with left
% p1(2) + alpha*dir(1) = 1
% int3 = p1 + alpha*dir
alpha = (1 - p1(1)) / dir(1);
intp{3} = p1 + alpha*dir;
intp{3}(1) = 1; % numerical precision issues....

% intersection with right
% p1(2) + alpha*dir(2) = imgwidth
% int4 = p1 + alpha*dir
alpha = (imgwidth - p1(1)) / dir(1);
intp{4} = p1 + alpha*dir;
intp{4}(1) = imgwidth; % numerical precision issues....

b(1) = isinimage(intp{1}, imgwidth, imgheight);
b(2) = isinimage(intp{2}, imgwidth, imgheight);
b(3) = isinimage(intp{3}, imgwidth, imgheight);
b(4) = isinimage(intp{4}, imgwidth, imgheight);

% should touch at least 2 edges
% greater than 2 when intersects exactly at corner of image
assert(b(1)+b(2)+b(3)+b(4) >= 2);  

b = find(b);
p1ext = intp{b(1)};
p2ext = intp{b(2)};


%%
function flag = isinimage(p, imgwidth, imgheight)
if p(1)>=1 && p(1)<=imgwidth && p(2)>=1 && p(2)<=imgheight
    flag = 1;
else
    flag = 0;
end

