function [p1ext p2ext degen] = extline(p1, p2, imgwidth, imgheight)
% p1: [x y]
% p2: [x y]

dir = p2 - p1;

% intersection with top
% p1(1) + alpha*dir(2) = 1
% int1 = p1 + alpha*dir
alpha = (1 - p1(2)) / dir(2);
intp{1} = p1 + alpha*dir;
intp{1}(2) = 1; % numerical precision issues....
intp{1} = intp{1}(:)'; % [2x1]

% intersection with bottom
% p1(1) + alpha*dir(2) = imgheight
% int2 = p1 + alpha*dir
alpha = (imgheight - p1(2)) / dir(2);
intp{2} = p1 + alpha*dir;
intp{2}(2) = imgheight; % numerical precision issues....
intp{2} = intp{2}(:)'; % [2x1]

% intersection with left
% p1(2) + alpha*dir(1) = 1
% int3 = p1 + alpha*dir
alpha = (1 - p1(1)) / dir(1);
intp{3} = p1 + alpha*dir;
intp{3}(1) = 1; % numerical precision issues....
intp{3} = intp{3}(:)'; % [2x1]

% intersection with right
% p1(2) + alpha*dir(2) = imgwidth
% int4 = p1 + alpha*dir
alpha = (imgwidth - p1(1)) / dir(1);
intp{4} = p1 + alpha*dir;
intp{4}(1) = imgwidth; % numerical precision issues....
intp{4} = intp{4}(:)'; % [2x1]

b(1) = isinimage(intp{1}, imgwidth, imgheight);
b(2) = isinimage(intp{2}, imgwidth, imgheight);
b(3) = isinimage(intp{3}, imgwidth, imgheight);
b(4) = isinimage(intp{4}, imgwidth, imgheight);

% should touch at least 2 edges
% greater than 2 when intersects exactly at corner of image
% assert(b(1)+b(2)+b(3)+b(4) >= 2);
if b(1)+b(2)+b(3)+b(4) < 2
    degen = 1;
    p1ext = [0 0];
    p2ext = [0 0];
    return;
else
    degen = 0;
end

% b = find(b);
% p1ext = intp{b(1)};
% p2ext = intp{b(2)};

p = cat(1, intp{logical(b)});
% p = unique(p, 'rows');
% assert(size(p,1)==2);
p1ext = p(1,:);
p2ext = p(2,:);


%%
function flag = isinimage(p, imgwidth, imgheight)
if p(1)>=1 && p(1)<=imgwidth && p(2)>=1 && p(2)<=imgheight
    flag = 1;
else
    flag = 0;
end

