% simple demo to convert a polygon mesh into a voxel representation

% code is 99.99% based on
% [1] http://www.mathworks.com/matlabcentral/fileexchange/24086-polygon2voxel
% [2] http://www.mathworks.com/matlabcentral/fileexchange/21044-3d-voxelizer

%   % Compile the c-coded function
%   mex polygon2voxel_double.c -v
clear all
close all
clc
load model;

vertices = vertices - repmat(mean(vertices,1),size(vertices,1),1);

FV.faces = faces;
s = 100/(max(max(vertices))-min(min(vertices)))
FV.vertices = (vertices-min(min(vertices)))*s;

Volume=polygon2voxel(FV,[100 100 100],'none');

%% visualization 1
figure
[X,Y,Z]=ind2sub(size(Volume),find(Volume(:)));
plot3(X,Y,Z,'.');
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
mean(X)
mean(Y)
mean(Z)

% %% visualization 2
% % 3d pirnter style visualization to add layer by layer
% plot3D(Volume,'timed', 0.1)
% 
% figure
% %% visualization 3
% for i=1:size(Volume,1)
%     imagesc(squeeze(Volume(i,:,:)));
%     axis equal;
%     axis tight;
%     axis off
%     title(i);
%     pause(0.1);
% end
% 
% for i=1:size(Volume,2)
%     imagesc(squeeze(Volume(:,i,:)));
%     axis equal;
%     axis tight;
%     axis off
%     title(i);
%     pause(0.1);
% end
% 
% for i=1:size(Volume,3)
%     imagesc(squeeze(Volume(:,:,i)));
%     axis equal;
%     axis tight;
%     axis off
%     title(i);
%     pause(0.1);
% end
