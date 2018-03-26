function quad = fit_bounding_quad(ptlist, ori, vp)


[pleft pright pbot2 ptop2 pbot3 ptop3 pleft3 pright3] = ...
    getextremepoints(ptlist, vp);

if ori==1
%     quad.pt(1,1:2) = line_intersect(ptop2, vp{2}, ptop3, vp{3}); % top
%     quad.pt(2,1:2) = line_intersect(ptop2, vp{2}, pbot3, vp{3}); % t2b3
%     quad.pt(3,1:2) = line_intersect(pbot2, vp{2}, ptop3, vp{3}); % b2t3
%     quad.pt(4,1:2) = line_intersect(pbot2, vp{2}, pbot3, vp{3}); % bot
    quad.pt(1,1:2) = line_intersect(ptop2, vp{2}, pleft3, vp{3}); % topleft
    quad.pt(2,1:2) = line_intersect(ptop2, vp{2}, pright3, vp{3}); % topright
    quad.pt(3,1:2) = line_intersect(pbot2, vp{2}, pleft3, vp{3}); % botleft
    quad.pt(4,1:2) = line_intersect(pbot2, vp{2}, pright3, vp{3}); % botleft
    if quad.pt(1,1)>quad.pt(2,1) % || quad.pt(3,1)>quad.pt(4,1)
        quad.pt = quad.pt([2 1 4 3],:);
    end
elseif ori==2
    quad.pt(1,1:2) = line_intersect(pleft,  vp{1}, ptop3, vp{3});  % topleft
    quad.pt(2,1:2) = line_intersect(pright, vp{1}, ptop3, vp{3}); % topright
    quad.pt(3,1:2) = line_intersect(pleft,  vp{1}, pbot3, vp{3});  % botleft
    quad.pt(4,1:2) = line_intersect(pright, vp{1}, pbot3, vp{3}); % botright
    if quad.pt(1,2)>quad.pt(3,2) % || quad.pt(2,2)>quad.pt(4,2)
        quad.pt = quad.pt([3 4 1 2],:);
    end
elseif ori==3
    quad.pt(1,1:2) = line_intersect(pleft,  vp{1}, ptop2, vp{2});  % topleft
    quad.pt(2,1:2) = line_intersect(pright, vp{1}, ptop2, vp{2}); % topright
    quad.pt(3,1:2) = line_intersect(pleft,  vp{1}, pbot2, vp{2});  % botleft
    quad.pt(4,1:2) = line_intersect(pright, vp{1}, pbot2, vp{2}); % botright
end


% % % % % % 
% % % % % pt = line_intersect(pleft,vp{1},ptop2,vp{2});
% % % % % junc2(1).pt = pt;
% % % % % junc2(1).type = 1;
% % % % % junc2(1).hclass = 2;
% % % % % junc2(1).vclass = 1;
% % % % % junc2(1).regionid = getregionid(pt(1),pt(2),vp);
% % % % % 
% % % % % pt = line_intersect(pright,vp{1},ptop2,vp{2});
% % % % % junc2(2).pt = pt;
% % % % % junc2(2).type = 2;
% % % % % junc2(2).hclass = 2;
% % % % % junc2(2).vclass = 1;
% % % % % junc2(2).regionid = getregionid(pt(1),pt(2),vp);
% % % % % 
% % % % % pt = line_intersect(pleft,vp{1},ptop3,vp{3});
% % % % % junc2(3).pt = pt;
% % % % % junc2(3).type = 1;
% % % % % junc2(3).hclass = 3;
% % % % % junc2(3).vclass = 1;
% % % % % junc2(3).regionid = getregionid(pt(1),pt(2),vp);
% % % % % 
% % % % % pt = line_intersect(pright,vp{1},ptop3,vp{3});
% % % % % junc2(4).pt = pt;
% % % % % junc2(4).type = 2;
% % % % % junc2(4).hclass = 3;
% % % % % junc2(4).vclass = 1;
% % % % % junc2(4).regionid = getregionid(pt(1),pt(2),vp);
% % % % % 
% % % % % %
% % % % % pt = line_intersect(pleft,vp{1},pbot2,vp{2});
% % % % % junc2(5).pt = pt;
% % % % % junc2(5).type = 3;
% % % % % junc2(5).hclass = 2;
% % % % % junc2(5).vclass = 1;
% % % % % junc2(5).regionid = getregionid(pt(1),pt(2),vp);
% % % % % 
% % % % % pt = line_intersect(pright,vp{1},pbot2,vp{2});
% % % % % junc2(6).pt = pt;
% % % % % junc2(6).type = 4;
% % % % % junc2(6).hclass = 2;
% % % % % junc2(6).vclass = 1;
% % % % % junc2(6).regionid = getregionid(pt(1),pt(2),vp);
% % % % % 
% % % % % pt = line_intersect(pleft,vp{1},pbot3,vp{3});
% % % % % junc2(7).pt = pt;
% % % % % junc2(7).type = 3;
% % % % % junc2(7).hclass = 3;
% % % % % junc2(7).vclass = 1;
% % % % % junc2(7).regionid = getregionid(pt(1),pt(2),vp);
% % % % % 
% % % % % pt = line_intersect(pright,vp{1},pbot3,vp{3});
% % % % % junc2(8).pt = pt;
% % % % % junc2(8).type = 4;
% % % % % junc2(8).hclass = 3;
% % % % % junc2(8).vclass = 1;
% % % % % junc2(8).regionid = getregionid(pt(1),pt(2),vp);
% % % % % 
% % % % % %
% % % % % pt = line_intersect(ptop2,vp{2},ptop3,vp{3});
% % % % % junc2(9).pt = pt;
% % % % % junc2(9).type = 1;
% % % % % junc2(9).hclass = 3;
% % % % % junc2(9).vclass = 2;
% % % % % junc2(9).regionid = getregionid(pt(1),pt(2),vp);
% % % % % 
% % % % % pt = line_intersect(pbot2,vp{2},pbot3,vp{3});
% % % % % junc2(10).pt = pt;
% % % % % junc2(10).type = 3;
% % % % % junc2(10).hclass = 3;
% % % % % junc2(10).vclass = 2;
% % % % % junc2(10).regionid = getregionid(pt(1),pt(2),vp);
% % % % % 
% % % % % % ***** TODO: junc11,12: type???
% % % % % pt = line_intersect(ptop2,vp{2},ptop3,vp{3});
% % % % % junc2(11).pt = pt;
% % % % % junc2(11).type = 2;
% % % % % junc2(11).hclass = 3;
% % % % % junc2(11).vclass = 2;
% % % % % junc2(11).regionid = getregionid(pt(1),pt(2),vp);
% % % % % 
% % % % % pt = line_intersect(pbot2,vp{2},pbot3,vp{3});
% % % % % junc2(12).pt = pt;
% % % % % junc2(12).type = 4;
% % % % % junc2(12).hclass = 3;
% % % % % junc2(12).vclass = 2;
% % % % % junc2(12).regionid = getregionid(pt(1),pt(2),vp);
% % % % % 
