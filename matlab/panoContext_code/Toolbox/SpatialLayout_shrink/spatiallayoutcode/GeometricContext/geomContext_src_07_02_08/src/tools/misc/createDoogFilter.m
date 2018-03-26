function DOOG = createDoogFilter(sigma,r,theta,gSize)
% DOOG = createDoogFilter(sigma,r,theta,gSize)
%
% Creates an difference of oriented gaussian filter (code from an old hw)
% sigma - lower sigmas result in "sharper" filters
% r - higher rs result in filters that are more highly directed
% gSize - the size of the filter to output
% theta - the angle of the direction in degrees
% Returns a DOOG filter

sigmaX=r*sigma;
sigmaY=sigma;
y0=sigma;
a=-1;
b=2;
c=1;

% create three gaussians
G1=createGaussian(0,y0,sigmaX,sigmaY,gSize);
G2=createGaussian(0,0,sigmaX,sigmaY,4*gSize+1);
G3=createGaussian(0,-y0,sigmaX,sigmaY,gSize);

% rotate the center gaussian according to theta
G2=imrotate(G2,theta,'bicubic','crop');
i1 = (4*gSize+1)/2 - gSize/2 + 1;
i2 = i1 + gSize - 1;
G2=G2(i1:i2,i1:i2);
G2=G2./sum(sum(G2));

% the final filter is a sum of the three
DOOG=a*G1+b*G2+c*G3;
DOOG=DOOG./sum(sum(DOOG));