function [ lines_ali ] = refitLineSegmentB( lines, vp, vpweight )
%REFITLINESEGMENTA refit direction of line segments 
%   lines: original line segments
%   vp: vannishing point
%   vpweight: if set to 0, lines will not change; if set to inf, lines will
%   be forced to pass vp
if ~exist('vpweight', 'var')
    vpweight = 0.1;
end

numSample = 100;
numLine = size(lines,1);
xyz = zeros(numSample+1,3);
wei = ones(numSample+1,1); wei(numSample+1) = vpweight*numSample;
lines_ali = lines;
for i = 1:numLine
    n = lines(i,1:3);
    sid = lines(i,5)*2*pi;
    eid = lines(i,6)*2*pi;
    if eid<sid
        x = linspace(sid,eid+2*pi,numSample);
        x = rem(x,2*pi);
%         x = sid-1:(eid-1+numBins);
%         x = rem(x,numBins) + 1;
    else
        x = linspace(sid,eid,numSample);
    end
%     u = -pi + (x'-1)*uBinSize + uBinSize/2; 
    u = -pi + x';
    v = computeUVN(n, u, lines(i,4));
    xyz(1:numSample,:) = uv2xyzN([u v], lines(i,4));
    xyz(numSample+1,:) = vp;
    [ ~, outputNM ] = curveFitting( xyz, wei );
    lines_ali(i,1:3) = outputNM;
end

end



