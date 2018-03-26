function [ rotImg, R ] = rotatePanorama( img, vp, R )
%ROTATEPANORAMA Rotate panorama
%   if R is given, vp (vanishing point) will be overlooked
%   otherwise R is computed from vp

[sphereH, sphereW, C] = size(img);
% rotImg = zeros( sphereH, sphereW, C);

%% new uv coordinates
[TX TY] = meshgrid(1:sphereW, 1:sphereH);
TX = TX(:);
TY = TY(:);
ANGx = (TX- sphereW/2 -0.5)/sphereW * pi *2 ;
ANGy = -(TY- sphereH/2 -0.5)/sphereH * pi;
uvNew = [ANGx ANGy];
xyzNew = uv2xyzN(uvNew,1);

%% rotation matrix
if nargin<3
    R = diag([1 1 1])/(vp');
end

xyzOld = (R\xyzNew')';
uvOld = xyz2uvN(xyzOld, 1);

% Px = uvOld(:,1)/2/pi*sphereW + 0.5 + sphereW/2;
% Py = -uvOld(:,2)/pi*sphereH + 0.5 + sphereH/2;
Px = (uvOld(:,1)+pi) / (2*pi) * sphereW + 0.5;
Py = (-uvOld(:,2) + pi/2) / pi * sphereH + 0.5;

Px = reshape(Px, [sphereH sphereW]);
Py = reshape(Py, [sphereH sphereW]);

% boundary
imgNew = double(zeros(sphereH+2, sphereW+2, C));
imgNew(2:end-1, 2:end-1, :) = img;
imgNew(2:end-1,1,:) = img(:,end,:);
imgNew(2:end-1,end,:) = img(:,1,:);
imgNew(1,2:sphereW/2+1,:) = img(1,sphereW:-1:sphereW/2+1,:);
imgNew(1,sphereW/2+2:end-1,:) = img(1,sphereW/2:-1:1,:);
imgNew(end,2:sphereW/2+1,:) = img(end,sphereW:-1:sphereW/2+1,:);
imgNew(end,sphereW/2+2:end-1,:) = img(1,sphereW/2:-1:1,:);
imgNew(1,1,:) = img(1,1,:);
imgNew(end,end,:) = img(end,end,:);
imgNew(1,end,:) = img(1,end,:);
imgNew(end,1,:) = img(end,1,:);

rotImg = warpImageFast(imgNew, Px+1, Py+1);
% rotImg = warpImageFast(img, Px, Py);

end

