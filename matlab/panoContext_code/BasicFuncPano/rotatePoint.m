function [ op ] = rotatePoint( p, R )
%ROTATEPOINT Rotate points
%   Detailed explanation goes here
op = (R * p')';
end

