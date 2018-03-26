function [ coor, tri ] = getUniformVector( level )
%GETUNIFORMVECTOR Summary of this function goes here
%   Detailed explanation goes here
if exist(sprintf('./BasicFuncPano/icosahedron2sphere/uniformvector_lvl%d.mat', level), 'file')
    load(sprintf('./BasicFuncPano/icosahedron2sphere/uniformvector_lvl%d.mat', level));
else
    [coor,tri] = icosahedron2sphere(level);
end

end

