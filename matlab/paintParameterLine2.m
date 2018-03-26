function [ panoEdgeC ] = paintParameterLine(parameterLine, width, height, img, c )
%PAINTPARAMETERLINE Paint parameterized line
%   Detailed explanation goes here
lines = parameterLine;

if nargin<=3
    panoEdgeC = zeros(height, width);
else
    img = im2double(img);
    panoEdgeC = imresize(img, [height width]);
    %if size(img,3)==3
    %    panoEdgeC = 0.5*rgb2gray(panoEdgeC);
    %end
end
% valid = true(size(lines,1),1);
% uv_vp = xyz2vp([vp;-vp]);
% vpm = min(floor( (uv_vp(:,1)-(-pi)) /(2*pi)*width)+1, width);
% vpn = min(floor( ((pi/2)-uv_vp(:,2))/(pi)*height )+1, height);
% valid = lines(:,4)==1;
% lines = lines(~valid,:); %horizontal
% lines = lines(valid,:); %vertical
num_sample = max(height,width);
for i = 1:size(lines,1)
%     fprintf('%d\n',i);
    n = lines(i,1:3);
    sid = lines(i,5)*2*pi;
    eid = lines(i,6)*2*pi;
    if eid<sid
        x = linspace(sid,eid+2*pi,num_sample);
        x = rem(x,2*pi);
%         x = sid-1:(eid-1+numBins);
%         x = rem(x,numBins) + 1;
    else
        x = linspace(sid,eid,num_sample);
    end
%     u = -pi + (x'-1)*uBinSize + uBinSize/2; 
    u = -pi + x';
    v = computeUVN(n, u, lines(i,4));
    xyz = uv2xyzN([u v], lines(i,4));
    uv = xyz2uvN( xyz, 1);
    
    m = min(floor( (uv(:,1)-(-pi)) /(2*pi)*width)+1, width);
    n = min(floor( ((pi/2)-uv(:,2))/(pi)*height )+1, height);
    %drawId = sub2ind([height width], n, m);
    %keyboard
    drawId = sub2ind([height width 3], n, m, repmat(c, size(m)));
    panoEdgeC(drawId) = i;
end

end

