function [ panoSegment ] = gbPanoSegment( img, sigma, k, minSz )
%GBPANOSEGMENT Summary of this function goes here
%   Detailed explanation goes here
% global config

[height, width, ~] = size(img);
img_smooth = smooth(img, sigma);

%% pixel segmentation, finish
% edges = zeros(width*height*4,3);
% [gridX, gridY] = meshgrid(1:width,1:height);
% num = 0;
% 
% x = gridX(:);
% y = gridY(:);
% 
% % xyz = uv2xyzN(coords2uv([x y], 2048, 1024));
% % angleNorm = 2*pi/width;
% 
% vector = [1 0; 0 1; 1 1; 1 -1];
% inda = sub2ind([height width], y, x);
% 
% for vid = 1:size(vector,1)
%     xv = x+vector(vid,1);
%     yv = y+vector(vid,2);
%     valid = xv>=1 & xv<=width & yv>=1 & yv<=height;
%     indav = inda(valid);
%     indbv = sub2ind([height width], yv(valid), xv(valid));
%     diff = (img_smooth(indav)-img_smooth(indbv)).^2 ...
%          + (img_smooth(indav+width*height)-img_smooth(indbv+width*height)).^2 ...
%          + (img_smooth(indav+2*width*height)-img_smooth(indbv+2*width*height)).^2;
% %     edges(num+1:num+sum(valid),:) = [x(valid)-1 + (y(valid)-1)*width ...
% %                                      xv(valid)-1 + (yv(valid)-1)*width ...
% %                                      sqrt(diff)];
% %     angleWeight = sqrt(1-dot(xyz(indav,:), xyz(indbv,:), 2).^2)/angleNorm;
%     angleWeight = ones(length(indav),1);
%     edges(num+1:num+sum(valid),:) = [indav indbv sqrt(diff).*angleWeight];
%     num = num + sum(valid);
% end
% 
% % connect left to right
% XY(1).xyLHS = [ones(height,1) (1:height)'];
% XY(1).xyRHS = [width*ones(height,1) (1:height)'];
% XY(2).xyLHS = [ones(height-1,1) (1:height-1)'];
% XY(2).xyRHS = [width*ones(height-1,1) (2:height)'];
% XY(3).xyLHS = [ones(height-1,1) (2:height)'];
% XY(3).xyRHS = [width*ones(height-1,1) (1:height-1)'];
% % 
% % combos = combntns(1:width,2);
% % XY(4).xyLHS = [combos(:,1) ones(size(combos,1),1)];
% % XY(4).xyRHS = [combos(:,2) ones(size(combos,1),1)];
% % XY(5).xyLHS = [combos(:,1) height*ones(size(combos,1),1)];
% % XY(5).xyRHS = [combos(:,2) height*ones(size(combos,1),1)];
% % 
% for sid = 1:length(XY)
%     xyLHS = XY(sid).xyLHS;
%     xyRHS = XY(sid).xyRHS;
%     indav = sub2ind([height width], xyLHS(:,2), xyLHS(:,1));
%     indbv = sub2ind([height width], xyRHS(:,2), xyRHS(:,1));
%     diff = (img_smooth(indav)-img_smooth(indbv)).^2 ...
%              + (img_smooth(indav+width*height)-img_smooth(indbv+width*height)).^2 ...
%              + (img_smooth(indav+2*width*height)-img_smooth(indbv+2*width*height)).^2;
%     addNum = size(xyLHS,1);
% %     edges(num+1:num+addNum,:) = [xyLHS(:,1)-1+(xyLHS(:,2)-1)*width ...
% %                                  xyRHS(:,1)-1+(xyRHS(:,2)-1)*width ...
% %                                  sqrt(diff)];
% %     angleWeight = sqrt(1-dot(xyz(indav,:), xyz(indbv,:), 2).^2)/angleNorm;
%     angleWeight = ones(length(indav),1);
%     edges(num+1:num+addNum,:) = [indav indbv sqrt(diff).*angleWeight];
%     num = num + addNum;
% end
% 
% edges = edges';
% segment = segmentGraphMex_edge(width*height, num, edges, k, minSz);
% L = unique(segment);
% temp = zeros(height, width);
% for i = 1:length(L)
%     temp(segment==L(i)) = i;
% end
% panoSegment = temp;
% 
% % segment = segmentGraphMex(width, height, num, edges, k, minSz);
% % L = unique(segment);
% % 
% % temp = zeros(size(segment));
% % for i = 1:length(L)
% %     temp(segment==L(i)) = i;
% % end

%% uniformly sample vectors on sphere and segment, test later
% global coor;
% global tri;
% load('./rectangleDetector/segmentation/uniformvector_lvl8.mat');
[coor, tri] = getUniformVector(8);
% load('./region_based_hypothesis/SketchTokens-master/models/forest/modelSmall.mat');
% st = stDetect( img, model );
% E = stToEdges( st, 1 );
[ E ] = getSketchTokenEdgemap( img );

% h = [0 -1 0; -1 4 -1; 0 -1 0];
% E = filter2(h,rgb2gray(img));
[EE, Ix, Iy] = dt2(double(E), 0.1, 0, 0.1, 0 );
% EE = zeros(1024, 2048);
% EE(747:846,1641:1697) = 1;

% [coor,tri] = icosahedron2sphere(level);
xySubs = uv2coords(xyz2uvN(coor), width, height);
xyinds = sub2ind([height width], xySubs(:,2), xySubs(:,1));
offset = width*height;

edges = [tri(:,1) tri(:,2); tri(:,2) tri(:,3); tri(:,3) tri(:,1)];
invert = edges(:,2)<edges(:,1);
edges(invert,:) = edges(invert,[2 1]);

uniEdges = unique(edges, 'rows');
weight = (img_smooth(xyinds(uniEdges(:,1)))-img_smooth(xyinds(uniEdges(:,2)))).^2 ...
       + (img_smooth(xyinds(uniEdges(:,1))+offset)-img_smooth(xyinds(uniEdges(:,2))+offset)).^2 ...
       + (img_smooth(xyinds(uniEdges(:,1))+2*offset)-img_smooth(xyinds(uniEdges(:,2))+2*offset)).^2;
gdweight = (EE(xyinds(uniEdges(:,1)))+EE(xyinds(uniEdges(:,2))))/2;
panoEdge = [uniEdges(:,1)'; uniEdges(:,2)'; sqrt(weight)'+10*double(gdweight)'];

maxID = size(coor,1);
num = size(uniEdges,1);

edgeLabel = segmentGraphMex_edge(maxID, num, panoEdge, k, minSz);

L = unique(edgeLabel);
temp = zeros(size(edgeLabel));
% center = zeros(length(L),2);
[gridX, gridY] = meshgrid(1:width, 1:height);
for i = 1:length(L)
%     fprintf('%d\n',i);
    temp(edgeLabel==L(i)) = i;
%     center(i,:) = [mean(gridX(edgeLabel==L(i))) mean(gridY(edgeLabel==L(i)))];
end


pixelvector = uv2xyzN(coords2uv([gridX(:) gridY(:)], width, height));

k = 1;
% [nnidx, dists] = annsearch( coor', pixelvector', k);
[nnidx, dists] = knnsearch( coor, pixelvector);

panoSegment = reshape(temp(nnidx), height, width);
% [ colors, pc ] = pcsel( center, 6 );
% C = [1 0 0;0 1 0;0 0 1;1 1 0;1 0 1;0 1 1];
% figure; imshow(label2rgb(panoSegment, C(colors,:)));

% panoSegment = zeros(height, width);
% for i = 1:length(L)
%     panoSegment(nnidx
% end


end

