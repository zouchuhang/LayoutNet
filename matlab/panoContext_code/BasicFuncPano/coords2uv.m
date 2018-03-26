function [ uv ] = coords2uv( coords, width, height )
%COORDS2UV Image coordinates (xy) to uv
%   Detailed explanation goes here
middleX = width/2+0.5;
middleY = height/2+0.5;
uv = [(coords(:,1)-middleX)./width*2*pi -(coords(:,2)-middleY)./height*pi];

end

