function [ out ] = computeUVN( n, in, planeID )
%COMPUTEUVN compute v given u and normal.
%   Detailed explanation goes here


if planeID==2
    n = [n(2) n(3) n(1)];
end
if planeID==3
    n = [n(3) n(1) n(2)];
end
bc = n(1)*sin(in) + n(2)*cos(in);
bs = n(3);
out = atan(-bc/bs);
end

