function warped_im = imgLookAt(im, CENTERx, CENTERy, new_imgH, fov )

%{
Citation:
J. Xiao, K. A. Ehinger, A. Oliva and A. Torralba.
Recognizing Scene Viewpoint using Panoramic Place Representation.
Proceedings of 25th IEEE Conference on Computer Vision and Pattern Recognition, 2012.
http://sun360.mit.edu
%}


sphereW = size(im,2);  sphereH = size(im,1);

warped_im = zeros(new_imgH,new_imgH,3);
[TXwarp TYwarp] = meshgrid(1:new_imgH, 1:new_imgH);
TX = TXwarp(:);                 TY = TYwarp(:);
TX = (TX -0.5 - new_imgH/2);    TY = (TY -0.5 - new_imgH/2);
% new_imgH = tan(fov/2) * R * 2
% TX = tan(ang/2) * R
r = new_imgH/2 / tan(fov/2);

% convert to 3D
R = sqrt(TY .^ 2 + r ^ 2);
ANGy = atan(- TY/r);
ANGy = ANGy + CENTERy;


X = sin(ANGy) .* R;
Y = - cos(ANGy) .* R;
Z = TX;

INDn = find(abs(ANGy) > pi/2); 

% project back to sphere

ANGx = atan(Z ./ -Y);
RZY = sqrt(Z .^2 + Y .^2);
ANGy = atan(X ./ RZY);
% INDn = find(abs(ANGy) > pi/2); 

ANGx(INDn) = ANGx(INDn)+pi; % if ANGy>pi/2, connect to +pi
ANGx = ANGx + CENTERx;

INDy = find(ANGy < -pi/2);  
ANGy(INDy) = - pi - ANGy(INDy) ;
ANGx(INDy) = ANGx(INDy) + pi;

INDx = find(ANGx <= -pi);     ANGx(INDx) = ANGx(INDx) + 2 * pi;
INDx = find(ANGx > pi);     ANGx(INDx) = ANGx(INDx) - 2 * pi;
INDx = find(ANGx > pi);     ANGx(INDx) = ANGx(INDx) - 2 * pi;
INDx = find(ANGx > pi);     ANGx(INDx) = ANGx(INDx) - 2 * pi;

% debug
% X: [-pi pi]
% Y: [-pi/2 pi/2]

Px = (ANGx+pi) / (2*pi) * sphereW + 0.5;
Py = ((- ANGy) + pi/2) / pi * sphereH + 0.5;

%Px(INDn)=1;
%Py(INDn)=1;

INDxx = find(Px<1);
Px(INDxx) = Px(INDxx) + sphereW;
im(:,sphereW+(1:2),:) = im(:,1:2,:);


% debug
% hold on
% plot(Px, Py, 'r.');


Px = reshape(Px, new_imgH, new_imgH);
Py = reshape(Py, new_imgH, new_imgH);

% finally warp image
warped_im = warpImageFast(im, Px, Py);

