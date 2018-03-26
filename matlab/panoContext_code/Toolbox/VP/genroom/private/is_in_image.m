function flag = is_in_image(p, imgwidth, imgheight, MARGIN)

% at least MARGIN pixels inside the image
if ~exist('MARGIN', 'var')
    MARGIN = 0;
end

if p(1)>=1+MARGIN && p(1)<=imgwidth-MARGIN && ...
   p(2)>=1+MARGIN && p(2)<=imgheight-MARGIN
    flag = 1;
else
    flag = 0;
end
