%demo for affine_fit
%Author: Adrien Leygue
%Date: December 3 2014
close all
clear all
figure;
%generate points that lie approximately in the Z=0 plane
N = 10;
[X,Y] = meshgrid(linspace(0,1,N));
XYZ_1 = [X(:) Y(:) 0.05*randn(N^2,1)];
plot3(XYZ_1(:,1),XYZ_1(:,2),XYZ_1(:,3),'r.');
hold on
%compute the normal to the plane and a point that belongs to the plane
[n_1,~,p_1] = affine_fit(XYZ_1);

%generate points that lie approximately in the Z=X plane
%the normal vector is
n_2_exact = [-sqrt(2)/2 0 sqrt(2)/2];
N = 12;
[X,Y] = meshgrid(linspace(0,1,N));
XYZ_2 = [X(:) Y(:) X(:)] + bsxfun(@times,0.05*randn(N^2,1),n_2_exact);
plot3(XYZ_2(:,1),XYZ_2(:,2),XYZ_2(:,3),'b.');


%compute the normal to the plane and a point that belongs to the plane
[n_2,V_2,p_2] = affine_fit(XYZ_2);

%plot the two points p_1 and p_2
plot3(p_1(1),p_1(2),p_1(3),'ro','markersize',15,'markerfacecolor','red');
plot3(p_2(1),p_2(2),p_2(3),'bo','markersize',15,'markerfacecolor','blue');

%plot the normal vector
quiver3(p_1(1),p_1(2),p_1(3),n_1(1)/3,n_1(2)/3,n_1(3)/3,'r','linewidth',2)
h = quiver3(p_2(1),p_2(2),p_2(3),n_2(1)/3,n_2(2)/3,n_2(3)/3,'b','linewidth',2)

%plot the two adjusted planes
[X,Y] = meshgrid(linspace(0,1,3));

%first plane
surf(X,Y, - (n_1(1)/n_1(3)*X+n_1(2)/n_1(3)*Y-dot(n_1,p_1)/n_1(3)),'facecolor','red','facealpha',0.5);

%second plane
%NB: if the plane is vertical the above method cannot be used, one should
%use the secont output of affine_fit which contains a base of the plane.
%this is illustrated below
%S1 and S2 are the coordinates of the plane points in the basis made of the
%columns ov V_2
[S1,S2] = meshgrid([-1 0 1]);
%generate the pont coordinates
X = p_2(1)+[S1(:) S2(:)]*V_2(1,:)';
Y = p_2(2)+[S1(:) S2(:)]*V_2(2,:)';
Z = p_2(3)+[S1(:) S2(:)]*V_2(3,:)';
%plot the plane
surf(reshape(X,3,3),reshape(Y,3,3),reshape(Z,3,3),'facecolor','blue','facealpha',0.5);

xlabel('x');
ylabel('y');
zlabel('z');
axis equal
%compute the angle between the planes in [0 90] degrees
angle = acosd(dot(n_1,n_2));
if angle>90
    angle = 180-angle;
end
angle

