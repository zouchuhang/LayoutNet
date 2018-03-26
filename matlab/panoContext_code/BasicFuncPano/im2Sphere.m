function [sphereImg validMap] = im2Sphere(im, imHoriFOV, sphereW, sphereH, x, y)
% perfect warping from separate image to panorama
% Work for x in [-pi,+pi], y in [-pi/2,+pi/2], and proper FOV
% For other (x,y,fov), it should also work, but depends on trigonometric 
% functions process in matlab, and not tested.

% map pixel in panorama to viewing direction
[TX TY] = meshgrid(1:sphereW, 1:sphereH);
TX = TX(:);
TY = TY(:);
ANGx = (TX- sphereW/2 -0.5)/sphereW * pi *2 ;
ANGy = -(TY- sphereH/2 -0.5)/sphereH * pi;

% compute the radius of ball
[imH, imW, ~] = size(im);
R = (imW/2) / tan(imHoriFOV/2);

% im is the tangent plane, contacting with ball at [x0 y0 z0]
x0 = R * cos(y) * sin(x);
y0 = R * cos(y) * cos(x);
z0 = R * sin(y);

% plane function: x0(x-x0)+y0(y-y0)+z0(z-z0)=0
% view line: x/alpha=y/belta=z/gamma
% alpha=cos(phi)sin(theta);  belta=cos(phi)cos(theta);  gamma=sin(phi)
alpha = cos(ANGy).*sin(ANGx);
belta = cos(ANGy).*cos(ANGx);
gamma = sin(ANGy);

% solve for intersection of plane and viewing line: [x1 y1 z1]
division = x0*alpha + y0*belta + z0*gamma;
x1 = R*R*alpha./division;
y1 = R*R*belta./division;
z1 = R*R*gamma./division;

% vector in plane: [x1-x0 y1-y0 z1-z0]
% positive x vector: vecposX = [cos(x) -sin(x) 0]
% positive y vector: vecposY = [x0 y0 z0] x vecposX
vec = [x1-x0 y1-y0 z1-z0];
vecposX = [cos(x) -sin(x) 0];
deltaX = (vecposX*vec') / sqrt(vecposX*vecposX');
vecposY = cross([x0 y0 z0], vecposX);
deltaY = (vecposY*vec') / sqrt(vecposY*vecposY');

% convert to im coordinates
Px = reshape(deltaX, [sphereH sphereW]) + (imW+1)/2;
Py = reshape(deltaY, [sphereH sphereW]) + (imH+1)/2;

% warp image
sphereImg = warpImageFast(im, Px, Py);
validMap = ~isnan(sphereImg(:,:,1));

% view direction: [alpha belta gamma]
% contacting point direction: [x0 y0 z0]
% so division>0 are valid region
validMap(division<0) = false;
validMap = repmat(validMap, [1 1 size(sphereImg,3)]);


