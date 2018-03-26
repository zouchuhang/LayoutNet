function [ lines_N ] = rotateLines( lines, R )
%ROTATELINES Rotate lines on panorama
%   Detailed explanation goes here

[numLine, dimLine] = size(lines);
lines_N = zeros(numLine, dimLine);
for i = 1:numLine
    n = lines(i,1:3);
    sid = lines(i,5)*2*pi-pi;
    eid = lines(i,6)*2*pi-pi;
    u = [sid;eid];
    v = computeUVN(n, u, lines(i,4));
    xyz = uv2xyzN([u v], lines(i,4));

    n_N = (R*n')'; n_N = n_N./norm(n_N);
    xyz_N = (R*xyz')';
    [~,lines_N(i,4)] = max(abs(n_N([3 1 2])));
    uv_N = xyz2uvN(xyz_N, lines_N(i,4));
    umax = max(uv_N(:,1))+pi;
    umin = min(uv_N(:,1))+pi;
    if umax-umin>pi
        lines_N(i,5:6) = [umax umin]/2/pi;
    else
        lines_N(i,5:6) = [umin umax]/2/pi;
    end
    
    lines_N(i,1:3) = n_N;   
%     lines_N(i,5:6) = (uv_N(:,1)'+pi)/2/pi;
    if dimLine>=7
        lines_N(i,7) = acos(dot(xyz_N(1,:), xyz_N(2,:), 2)./(norm(xyz_N(1,:))*norm(xyz_N(2,:))));
        % lines_N(i,7) = lines(i,7); % this should be ok as well
    end
    if dimLine>=8
        lines_N(i,8) = lines(i,8);
    end
    
end

end

