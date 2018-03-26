function [ b ] = judgeObbConnect( obb1, obb2, t )
%JUDGEOBBCONNECT not finished
%   Detailed explanation goes here

nBb1 = size(obb1,2);
nBb2 = size(obb2,2);

volume1 = cuboidVolume(obb1);
volume2 = cuboidVolume(obb2);
intersection = cuboidIntersectionVolume(obb1,obb2);
% union = repmat(volume1',1,nBb2)+repmat(volume2,nBb1,1)-intersection;
minvol = min(repmat(volume1',1,nBb2),repmat(volume2,nBb1,1));

scoreMatrix = intersection ./ minvol;
b = scoreMatrix>t;

end

