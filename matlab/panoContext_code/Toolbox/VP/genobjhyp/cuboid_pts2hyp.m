function cuboidhyp = cuboid_pts2hyp(pts, vp)

junc3dummy = struct('type',[],'pt',[],'regionid',[]);
cuboidhyp.junc3(8) = junc3dummy; % initialize

for i = 1:8
    cuboidhyp.junc3(i).type = i;
    cuboidhyp.junc3(i).pt = pts(i,:);
    cuboidhyp.junc3(i).regionid = getregionid(pts(i,1), pts(i,2), vp);
end
