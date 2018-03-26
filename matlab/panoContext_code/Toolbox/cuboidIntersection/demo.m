% written by jianxiong xiao
% the 3D cuboids should align with z axis
% each column is a 3D cuboids
% [x1 y1 x2 y2 x3 y3 x4 y4 zMin zMax]'
% this function will compute the intersection volume divided by the union volume
% between two groups of cuboids

% to compile to:
% mex gpc.c cuboidIntersectionVolume.c -O -output cuboidIntersectionVolume

bb1 = [ ...
         0    1.6303;
         0    2.0759;
    1.0000    1.5742;
         0    2.0611;
    1.0000    1.1589;
    1.0000    3.6412;
         0    1.2150;
    1.0000    3.6560;
         0   -1.1167;
    1.0000    1.2677];

bb2 = [ ...
    0.8911    1.2217    1.6303;
    1.8830    3.6303    2.0759;
   -1.2182   -1.4967    1.5742;
    1.3287    2.9274    2.0611;
   -1.6293   -1.5033    1.1589;
    2.8928    2.9531    3.6412;
    0.4800    1.2150    1.2150;
    3.4472    3.6560    3.6560;
    1.1863   -1.1667   -1.1167;
    1.2512    1.2677    1.2677];
    
    
nBb1 = size(bb1,2);
nBb2 = size(bb2,2);

volume1 = cuboidVolume(bb1);
volume2 = cuboidVolume(bb2);
intersection = cuboidIntersectionVolume(bb1,bb2);
union = repmat(volume1',1,nBb2)+repmat(volume2,nBb1,1)-intersection;

scoreMatrix = intersection ./ union;

scoreMatrix