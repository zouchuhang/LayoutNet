function [ out3D, out2D, valid, division ] = xyz2view( xyz, xyz0, imH, imW )
%XYZ2VIEW 3D points project to 2D image
%   xyz: 3D points
%   xyz0: 3D position of the center of the view
R = sqrt(sum(xyz0.^2));

% im is the tangent plane, contacting with ball at [x0 y0 z0]
x0 = xyz0(1);
y0 = xyz0(2);
z0 = xyz0(3);

uv = xyz2uvN(xyz, 1);
ANGx = uv(:,1); ANGy = uv(:,2);
uv = xyz2uvN(xyz0,1);
x = uv(1);  y = uv(2);
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
deltaX = (vecposX*vec') / sqrt(vecposX*vecposX') + (imW+1)/2;
vecposY = cross([x0 y0 z0], vecposX);
deltaY = (vecposY*vec') / sqrt(vecposY*vecposY') + (imH+1)/2;

out3D = [x1 y1 z1];
out2D = [deltaX' deltaY'];
valid = division>-0.000001;
end

