function [ xyz ] = uv2xyzN( uv, planeID )
%UV2XYZN Summary of this function goes here
%   Detailed explanation goes here
if ~exist('planeID','var')
    planeID = 1;
end

ID1 = rem(planeID-1+0,3)+1;
ID2 = rem(planeID-1+1,3)+1;
ID3 = rem(planeID-1+2,3)+1;

xyz = zeros(size(uv,1),3);
xyz(:,ID1) = cos(uv(:,2)).*sin(uv(:,1));
xyz(:,ID2) = cos(uv(:,2)).*cos(uv(:,1));
xyz(:,ID3) = sin(uv(:,2));

end

