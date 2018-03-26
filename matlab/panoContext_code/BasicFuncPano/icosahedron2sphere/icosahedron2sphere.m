function [coor,tri] = icosahedron2sphere(level)

% copyright by Jianxiong Xiao http://mit.edu/jxiao
% this function use a icosahedron to sample uniformly on a sphere
%{
Please cite this paper if you use this code in your publication:
J. Xiao, T. Fang, P. Zhao, M. Lhuillier, and L. Quan
Image-based Street-side City Modeling
ACM Transaction on Graphics (TOG), Volume 28, Number 5
Proceedings of ACM SIGGRAPH Asia 2009
%}


a= 2/(1+sqrt(5));
M=[
    0 a -1 a 1 0 -a 1 0
    0 a 1 -a 1 0 a 1 0
    0 a 1 0 -a 1 -1 0 a
    0 a 1 1 0 a 0 -a 1
    0 a -1 0 -a -1 1 0 -a
    0 a -1 -1 0 -a 0 -a -1
    0 -a 1 a -1 0 -a -1 0
    0 -a -1 -a -1 0 a -1 0
    -a 1 0 -1 0 a -1 0 -a
    -a -1 0 -1 0 -a -1 0 a
    a 1 0 1 0 -a 1 0 a
    a -1 0 1 0 a 1 0 -a
    0 a 1 -1 0 a -a 1 0
    0 a 1 a 1 0 1 0 a
    0 a -1 -a 1 0 -1 0 -a
    0 a -1 1 0 -a a 1 0
    0 -a -1 -1 0 -a -a -1 0
    0 -a -1 a -1 0 1 0 -a
    0 -a 1 -a -1 0 -1 0 a
    0 -a 1 1 0 a a -1 0
    ];

coor = reshape(M',3,60)';
%[M(:,[1 2 3]); M(:,[4 5 6]); M(:,[7 8 9])];


[coor, ~, idx] = unique(coor,'rows');

tri = reshape(idx,3,20)';

%{
for i=1:size(tri,1)
    x(1)=coor(tri(i,1),1);
    x(2)=coor(tri(i,2),1);
    x(3)=coor(tri(i,3),1);
    y(1)=coor(tri(i,1),2);
    y(2)=coor(tri(i,2),2);
    y(3)=coor(tri(i,3),2);
    z(1)=coor(tri(i,1),3);
    z(2)=coor(tri(i,2),3);
    z(3)=coor(tri(i,3),3);
    patch(x,y,z,'r');
end

axis equal
axis tight
%}

% extrude
coor = coor ./ repmat(sqrt(sum(coor .* coor,2)),1, 3);

for i=1:level
    m = 0;
    for t=1:size(tri,1)
        n = size(coor,1);
        coor(n+1,:) = ( coor(tri(t,1),:) + coor(tri(t,2),:) ) / 2;
        coor(n+2,:) = ( coor(tri(t,2),:) + coor(tri(t,3),:) ) / 2;
        coor(n+3,:) = ( coor(tri(t,3),:) + coor(tri(t,1),:) ) / 2;
        
        triN(m+1,:) = [n+1     tri(t,1)    n+3];
        triN(m+2,:) = [n+1     tri(t,2)    n+2];
        triN(m+3,:) = [n+2     tri(t,3)    n+3];
        triN(m+4,:) = [n+1     n+2         n+3];
        
        n = n+3;
        
        m = m+4;
        
    end
    tri = triN;
    
    % uniquefy
    [coor, ~, idx] = unique(coor,'rows');
    tri = idx(tri);
    
    % extrude
    coor = coor ./ repmat(sqrt(sum(coor .* coor,2)),1, 3);
end

% vertex number: 12  42  162  642
