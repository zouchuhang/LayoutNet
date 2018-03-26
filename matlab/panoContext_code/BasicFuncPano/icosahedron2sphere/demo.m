% copyright by Jianxiong Xiao http://mit.edu/jxiao
% demo how to sample uniformly on a sphere

%{
Please cite this paper if you use this code in your publication:
J. Xiao, T. Fang, P. Zhao, M. Lhuillier, and L. Quan
Image-based Street-side City Modeling
ACM Transaction on Graphics (TOG), Volume 28, Number 5
Proceedings of ACM SIGGRAPH Asia 2009
%}

clear
clc

tic;
points = icosahedron2sphere(0);
toc;
subplot(1,4,1);
plot3(points(:,1),points(:,2),points(:,3),'.')
title(sprintf('Level %d with %d points',0,size(points,1)))
axis equal
axis tight

tic;
points = icosahedron2sphere(1);
toc;
subplot(1,4,2);
plot3(points(:,1),points(:,2),points(:,3),'.')
title(sprintf('Level %d with %d points',1,size(points,1)))
axis equal
axis tight

tic;
points = icosahedron2sphere(2);
toc;
subplot(1,4,3);
plot3(points(:,1),points(:,2),points(:,3),'.')
title(sprintf('Level %d with %d points',2,size(points,1)))
axis equal
axis tight

tic;
points = icosahedron2sphere(4);
toc;
subplot(1,4,4);
plot3(points(:,1),points(:,2),points(:,3),'.')
title(sprintf('Level %d with %d points',3,size(points,1)))
axis equal
axis tight
