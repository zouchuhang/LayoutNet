function [ lines_vp ] = getVanishingLine( lines, vp, d_tol )
%GETVANISHINGLINE filter out line segments with vannishing points
%   line segments in lines will be saved if only within d_tol degree with
%   vp
nml = lines(:,1:3);
cosint = dot(nml, repmat(vp, [size(nml,1) 1]), 2);
valid = abs(cosint)<cos((90-d_tol)*pi/180);
lines_vp = lines(valid,:);

numLine = size(lines_vp,1);
valid = true(numLine,1);
for i = 1:numLine
    us = lines_vp(i,5);
    ue = lines_vp(i,6);
    u = [us;ue]*2*pi-pi;
    v = computeUVN(lines_vp(i,1:3), u, lines_vp(i,4));
    xyz = uv2xyzN([u v], lines_vp(i,4));
    x = linspace(xyz(1,1),xyz(2,1),100);
    y = linspace(xyz(1,2),xyz(2,2),100);
    z = linspace(xyz(1,3),xyz(2,3),100);
    xyz = [x' y' z'];
    xyz = xyz ./ repmat(sqrt(sum(xyz.^2,2)),[1 3]);
    ang = acos( abs(dot(xyz, repmat(vp, [100 1]), 2)));
    valid(i) = ~any(ang<10*pi/180);
end
lines_vp = lines_vp(valid,:);

end

