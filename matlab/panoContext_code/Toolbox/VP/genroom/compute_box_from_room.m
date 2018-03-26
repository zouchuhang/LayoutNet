function roomhyp = compute_box_from_room(roomhyp, vp, imgwidth, imgheight)
% box(1): left wall
% box(2): middle wall
% box(3): right wall
% box(4): ceiling
% box(5): floor

% image boundary polygon
imgpoly = [1 1; 1 imgheight; imgwidth imgheight; imgwidth 1; 1 1];

for i = 1:length(roomhyp)
    rh = roomhyp(i);

    if rh.type==1
        pt = roomhyp(i).corner(3).pt;
        boxpoly(1).p1 = [];
        boxpoly(1).p2 = extend_to_imageedge(pt, vp{3}, 'away', imgwidth, imgheight);
        boxpoly(1).p3 = extend_to_imageedge(pt, vp{1}, 'up', imgwidth, imgheight);
        boxpoly(1).p4 = pt;
        boxpoly(1).orient = 2;
        
        pt = roomhyp(i).corner(3).pt;
        boxpoly(2).p1 = extend_to_imageedge(pt, vp{1}, 'up', imgwidth, imgheight);
        boxpoly(2).p2 = pt;
        boxpoly(2).p3 = [];
        boxpoly(2).p4 = extend_to_imageedge(pt, vp{2}, 'right', imgwidth, imgheight);
        boxpoly(2).orient = 3;
        
    elseif rh.type==2
        pt = roomhyp(i).corner(4).pt;
        boxpoly(1).p1 = [];
        boxpoly(1).p2 = extend_to_imageedge(pt, vp{2}, 'left', imgwidth, imgheight);
        boxpoly(1).p3 = extend_to_imageedge(pt, vp{1}, 'up', imgwidth, imgheight);
        boxpoly(1).p4 = pt;
        boxpoly(1).orient = 3;
        
        pt = roomhyp(i).corner(4).pt;
        boxpoly(2).p1 = extend_to_imageedge(pt, vp{1}, 'up', imgwidth, imgheight);
        boxpoly(2).p2 = pt;
        boxpoly(2).p3 = [];
        boxpoly(2).p4 = extend_to_imageedge(pt, vp{3}, 'away', imgwidth, imgheight);
        boxpoly(2).orient = 2;
    elseif rh.type==3
        boxpoly(1).p1 = extend_to_imageedge(roomhyp(i).corner(1).pt, vp{3}, 'away', imgwidth, imgheight);
        boxpoly(1).p2 = extend_to_imageedge(roomhyp(i).corner(3).pt, vp{3}, 'away', imgwidth, imgheight);
        boxpoly(1).p3 = roomhyp(i).corner(1).pt;
        boxpoly(1).p4 = roomhyp(i).corner(3).pt;
        boxpoly(1).orient = 2;
        
        boxpoly(2).p1 = roomhyp(i).corner(1).pt;
        boxpoly(2).p2 = roomhyp(i).corner(3).pt;
        boxpoly(2).p3 = extend_to_imageedge(roomhyp(i).corner(1).pt, vp{2}, 'right', imgwidth, imgheight);
        boxpoly(2).p4 = extend_to_imageedge(roomhyp(i).corner(3).pt, vp{2}, 'right', imgwidth, imgheight);
        boxpoly(2).orient = 3;
    elseif rh.type==4
        boxpoly(1).p1 = extend_to_imageedge(roomhyp(i).corner(2).pt, vp{2}, 'left', imgwidth, imgheight);
        boxpoly(1).p2 = extend_to_imageedge(roomhyp(i).corner(4).pt, vp{2}, 'left', imgwidth, imgheight);
        boxpoly(1).p3 = roomhyp(i).corner(2).pt;
        boxpoly(1).p4 = roomhyp(i).corner(4).pt;
        boxpoly(1).orient = 3;
        
        boxpoly(2).p1 = roomhyp(i).corner(2).pt;
        boxpoly(2).p2 = roomhyp(i).corner(4).pt;
        boxpoly(2).p3 = extend_to_imageedge(roomhyp(i).corner(2).pt, vp{3}, 'away', imgwidth, imgheight);
        boxpoly(2).p4 = extend_to_imageedge(roomhyp(i).corner(4).pt, vp{3}, 'away', imgwidth, imgheight);
        boxpoly(2).orient = 2;
    elseif rh.type==5
        boxpoly(1).p1 = extend_to_imageedge(roomhyp(i).corner(1).pt, vp{3}, 'away', imgwidth, imgheight);
        boxpoly(1).p2 = extend_to_imageedge(roomhyp(i).corner(3).pt, vp{3}, 'away', imgwidth, imgheight);
        boxpoly(1).p3 = roomhyp(i).corner(1).pt;
        boxpoly(1).p4 = roomhyp(i).corner(3).pt;
        boxpoly(1).orient = 2;
        
        boxpoly(2).p1 = roomhyp(i).corner(1).pt;
        boxpoly(2).p2 = roomhyp(i).corner(3).pt;
        boxpoly(2).p3 = roomhyp(i).corner(2).pt;
        boxpoly(2).p4 = roomhyp(i).corner(4).pt;
        boxpoly(2).orient = 3;
        
        boxpoly(3).p1 = roomhyp(i).corner(2).pt;
        boxpoly(3).p2 = roomhyp(i).corner(4).pt;
        boxpoly(3).p3 = extend_to_imageedge(roomhyp(i).corner(2).pt, vp{3}, 'away', imgwidth, imgheight);
        boxpoly(3).p4 = extend_to_imageedge(roomhyp(i).corner(4).pt, vp{3}, 'away', imgwidth, imgheight);
        boxpoly(3).orient = 2;
        
    elseif rh.type==6
        boxpoly(1).p1 = extend_to_imageedge(roomhyp(i).corner(1).pt, vp{3}, 'away', imgwidth, imgheight);
        boxpoly(1).p2 = [];
        boxpoly(1).p3 = roomhyp(i).corner(1).pt;
        boxpoly(1).p4 = extend_to_imageedge(roomhyp(i).corner(1).pt, vp{1}, 'down', imgwidth, imgheight);
        boxpoly(1).orient = 2;
        
        boxpoly(2).p1 = roomhyp(i).corner(1).pt;
        boxpoly(2).p2 = extend_to_imageedge(roomhyp(i).corner(1).pt, vp{1}, 'down', imgwidth, imgheight);
        boxpoly(2).p3 = extend_to_imageedge(roomhyp(i).corner(1).pt, vp{2}, 'right', imgwidth, imgheight);
        boxpoly(2).p4 = [];
        boxpoly(2).orient = 3;
    elseif rh.type==7
        boxpoly(1).p1 = extend_to_imageedge(roomhyp(i).corner(2).pt, vp{2}, 'left', imgwidth, imgheight);
        boxpoly(1).p2 = [];
        boxpoly(1).p3 = roomhyp(i).corner(2).pt;
        boxpoly(1).p4 = extend_to_imageedge(roomhyp(i).corner(2).pt, vp{1}, 'down', imgwidth, imgheight);
        boxpoly(1).orient = 3;
        
        boxpoly(2).p1 = roomhyp(i).corner(2).pt;
        boxpoly(2).p2 = extend_to_imageedge(roomhyp(i).corner(2).pt, vp{1}, 'down', imgwidth, imgheight);
        boxpoly(2).p3 = extend_to_imageedge(roomhyp(i).corner(2).pt, vp{3}, 'away', imgwidth, imgheight);
        boxpoly(2).p4 = [];
        boxpoly(2).orient = 2;
    end
    
	roomhyp(i).box = boxpoly;
end

%%
function ptedge = extend_to_imageedge(pt, vp, whichside, imgwidth, imgheight)
% whichside: 'up','down'(for vp{1}), 'left','right'(for vp{2}), 'away'(for vp{3})
pt = pt(:)';
vp = vp(:)';

dir = pt - vp;
if     strcmp(whichside, 'up'),    if sign(dir(2))==1,  dir = -dir; end
elseif strcmp(whichside, 'down'),  if sign(dir(2))==-1, dir = -dir; end
elseif strcmp(whichside, 'left'),  if sign(dir(1))==1,  dir = -dir; end
elseif strcmp(whichside, 'right'), if sign(dir(1))==-1, dir = -dir; end
% elseif strcmp(whichside, 'away'), % nothing
end
% if sign(dir(1))~=whichside
% 	dir = -dir;
% end

% intersection with top
% pt(2) + alpha*dir(2) = 1
% int1 = pt + alpha*dir
alpha(1) = (1 - pt(2)) / dir(2);
intp{1} = pt + alpha(1)*dir;
intp{1}(2) = 1; % numerical precision issues....

% intersection with bottom
% pt(2) + alpha*dir(2) = imgheight
% int2 = pt + alpha*dir
alpha(2) = (imgheight - pt(2)) / dir(2);
intp{2} = pt + alpha(2)*dir;
intp{2}(2) = imgheight; % numerical precision issues....

% intersection with left
% pt(1) + alpha*dir(1) = 1
% int3 = pt + alpha*dir
alpha(3) = (1 - pt(1)) / dir(1);
intp{3} = pt + alpha(3)*dir;
intp{3}(1) = 1; % numerical precision issues....

% intersection with right
% p1(1) + alpha*dir(1) = imgwidth
% int4 = pt + alpha*dir
alpha(4) = (imgwidth - pt(1)) / dir(1);
intp{4} = pt + alpha(4)*dir;
intp{4}(1) = imgwidth; % numerical precision issues....

b(1) = alpha(1)>0 && is_in_image(intp{1}, imgwidth, imgheight);
b(2) = alpha(2)>0 && is_in_image(intp{2}, imgwidth, imgheight);
b(3) = alpha(3)>0 && is_in_image(intp{3}, imgwidth, imgheight);
b(4) = alpha(4)>0 && is_in_image(intp{4}, imgwidth, imgheight);

b = find(b);
% assert(length(b)==1);
if length(b)==0
    ptedge = [];
elseif length(b)==1
    ptedge = intp{b};
else
    error('aosjfoiajeoijfoeaijfoewjfe');
end

