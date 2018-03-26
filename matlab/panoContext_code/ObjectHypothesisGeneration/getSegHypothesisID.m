function [ cuboid ] = getSegHypothesisID( labelMap, minSize, vanishing, ID )
%GETBINARYMAP Summary of this function goes here
%   Detailed explanation goes here
segment = labelMap;% + 1;
numSeg = round(max(segment(:)));
[H,W] = size(segment);

invalidMask = false(H,W);
invalidMask(:,1:10) = true;   invalidMask(:,W-9:W) = true;
invalidMask(1:10,:) = true;   invalidMask(H-9:H,:) = true;
   
fprintf('*1');
 
if nargin<=3
    n = histc(segment(:),1:numSeg);
    segValid = n>minSize;
else
    segValid = false(numSeg,1);
    segValid(ID) = true;
end
counter = 0;
% maxHypot= struct([]);
% hypot = zeros(6,2,1000);
% score = zeros(1,1000);

% binaryMaps = zeros(H,W,100);
SE = strel('disk', 5, 4);
SS = strel('disk', 2, 4);

fprintf('*2');

for i = find(segValid)'
    seg = segment == i;
    
    binaryMap = imdilate( seg, SE);
    binaryMap = imerode( binaryMap, SE);
    binaryMap = imerode( binaryMap, SS);
    binaryMap = imdilate( binaryMap, SS);
    binaryMap = imfill(binaryMap, 'hole');
    
    CC = bwconncomp(binaryMap);
    valid = true(CC.NumObjects,1);
    for j = 1:CC.NumObjects
        if length(CC.PixelIdxList{j})<minSize
            valid(j) = false;
        end
        if any(invalidMask(CC.PixelIdxList{j}))
            valid(j) = false;
        end
    end
    CC.NumObjects = sum(valid);
    CC.PixelIdxList = CC.PixelIdxList(valid);
    
%     R = regionprops( CC, 'BoundingBox', 'Area', 'ConvexArea', 'Centroid');
    for j = 1:CC.NumObjects
        subMask = false(H,W);
        subMask(CC.PixelIdxList{j}) = true;
        counter = counter + 1;
        maxHypot(counter) = region2cuboid( subMask, vanishing );

%         figure; imshow(subMask); hold on;
%         color = [1 0 0; 0.5 0 0; 0 1 0; 0 0.5 0; 0 0 1; 0 0 0.5];
%         points = maxHypot(counter).points;
%         for k = 1:6
%             line(points([1 2 3 4 5 6 1],1), points([1 2 3 4 5 6 1],2), 'Color', color(k,:), 'LineWidth', 2);
%         end
%         for k = 1:6
%             plot(points(k,1), points(k,2), '--rs','LineWidth',2,...
%                         'MarkerEdgeColor','k',...
%                         'MarkerFaceColor',color(k,:),...
%                         'MarkerSize',5);
%         end
       
    end
    
end

fprintf('*3');

if exist('maxHypot','var')
    score = [maxHypot.score];
    cuboid = maxHypot(score>0.75);
else
    cuboid = [];
end

fprintf('*4');

%%
% numCuboid = length(score);
% % rectangle = repmat(struct('xyzBox',zeros(1000,12),'xyBox',zeros(1000,8),'count',0),3,1);
% cuboid = struct('views', zeros(numCuboid,3), 'xyBox', zeros(numCuboid,24), 'score', zeros(numCuboid,1));
% for i = 1:numCuboid
%     points = hypot(:,:,i);
%     [x2,y2] = poly2cw(points(:,1), points(:,2));
%     points = [x2 y2];
%     
%     lines = points([2 3 4 5 6 1],:) - points([1 2 3 4 5 6],:);
%     direct = lines./repmat(sqrt(sum(lines.^2,2)), 1, 2);
%     midpoints = (points([2 3 4 5 6 1],:) + points([1 2 3 4 5 6],:))/2;
%     A = zeros(6,3);
%     for j = 1:3
%         vanline = [midpoints(:,1)-vanishing(j,1) midpoints(:,2)-vanishing(j,2)];
%         vandire = vanline./repmat(sqrt(sum(vanline.^2,2)), 1, 2);
%         A(:,j) = dot(direct, vandire, 2);
%     end
%     [~,direct_assign] = max( abs(A), [], 2);
%     % define type
%     type_indicate = direct_assign(1:3)==direct_assign(4:6);
%     if sum(type_indicate)==1 % two faces
%         I = find(type_indicate);
%         pointID1 = rem((I-1)+2,6)+1;
%         pointID2 = rem((pointID1-1)+3,6)+1;
%         testdir = points(pointID1,:)-points(pointID2,:);
%         testdir = testdir ./ norm(testdir);
%         refedir = vanishing(direct_assign(I),:) - (points(pointID1,:)+points(pointID2,:))/2;
%         refedir = refedir ./ norm(refedir);
%         if abs(dot(testdir, refedir))>0.95 % valid, connect line follow 3rd direction
%             pointID = rem([I I+1 I+2 I+5]-1,6)+1;
%             directIDs = direct_assign(pointID([4 1 2]));
%             directID = setdiff(1:3,directIDs);
% %             count1 = rectangle(directID).count + 1;
% %             rectangle(directID).xyBox(count1,:) = reshape(points(pointID,:)',1,[]);
% %             rectangle(directID).count = rectangle(directID).count + 1;
%             cuboid.xyBox(i,1:8) = reshape(points(pointID,:)',1,[]);
%             cuboid.views(i,1) = directID;
%             
%             pointID = rem([I+2 I+3 I+4 I+5]-1,6)+1;
%             directIDs = direct_assign(pointID([1 2 3]));
%             directID = setdiff(1:3,directIDs);
% %             count2 = rectangle(directID).count + 1;
% %             rectangle(directID).xyBox(count2,:) = reshape(points(pointID,:)',1,[]);
% %             rectangle(directID).count = rectangle(directID).count + 1;   
%             cuboid.xyBox(i,9:16) = reshape(points(pointID,:)',1,[]);
%             cuboid.views(i,2) = directID;
%             
% %             cuboid.combo(i,1:2) = [count1 count2];
%             cuboid.score(i) = score(i);
%         end
%     elseif sum(type_indicate)==3
%         % three faces
% %         [~,lowID] = max(points(:,2));
% %         ptrID = rem( lowID+[-2 0 +2 -1 +1 +3]-1, 6) + 1;
% %         ptr = zeros(7,2);
% %         ptr(2:7,:) = points(ptrID,:);
% %         [point2D, point3D, K] = fitCuboid(ptr(:,[2 1]), 1);
% %         points(ptrID,:) = point2D([2 1],2:7)';
% 
% %         pointIDs = zeros(1,3);
% %         directIDs = zeros(1,3);
%          I = find(direct_assign==3);
%          verticalID = rem([I(1)+2 I(1)+5]-1,6) + 1;
%         if norm(points(verticalID(1),:)-vanishing(3,:))<norm(points(verticalID(2),:)-vanishing(3,:))
%             pointIDs = rem([verticalID(1) verticalID(1)+2 verticalID(1)+4]-1,6)+1;
%         else
%             pointIDs = rem([verticalID(2) verticalID(2)+2 verticalID(2)+4]-1,6)+1;
%         end
%          
% %         if points(verticalID(1),2)<points(verticalID(2),2)
% %             pointIDs = rem([verticalID(2) verticalID(2)+2 verticalID(2)+4]-1,6)+1;
% %         else
% %             pointIDs = rem([verticalID(1) verticalID(1)+2 verticalID(1)+4]-1,6)+1;
% %         end
%         directIDs = direct_assign(rem(pointIDs-1+1,6)+1);
% %         for j = 1:3
% %             directID = direct_assign( j+1 );
% %             directIDs(j) = directID;
% %             if norm(points(j,:)-vanishing(directID,:))<norm(points(j+3,:)-vanishing(directID,:))
% %                 pointIDs(j) = j;
% %             else
% %                 pointIDs(j) = j+3;
% %             end
% %         end
%         lines = [points(pointIDs(1),:) vanishing(directIDs(1),:); ...
%                  points(pointIDs(2),:) vanishing(directIDs(2),:); ...
%                  points(pointIDs(3),:) vanishing(directIDs(3),:)];
%         intersect = intersectPoint(lines([1 2 3],:), lines([2 3 1],:));
% %         if sum(sum((intersect([1 2 3],:)-intersect([2 3 1],:)).^2,2))<10
%         
%             addPoint = sum(intersect,1)./3;
%             
% %             count1 = rectangle(directIDs(1)).count + 1;
% %             rectangle(directIDs(1)).count = rectangle(directIDs(1)).count + 1;
% %             rectangle(directIDs(1)).xyBox(count1,:) = reshape([points([1 2 3],:);addPoint]',1,[]);
%             cuboid.xyBox(i,1:8) = reshape([points(rem(pointIDs(1)+[0 1 2]-1,6)+1,:);addPoint]',1,[]);
%             cuboid.views(i,1) = setdiff(1:3, direct_assign([1 2]));
% %             count2 = rectangle(directIDs(2)).count + 1;
% %             rectangle(directIDs(2)).count = rectangle(directIDs(2)).count + 1;
% %             rectangle(directIDs(2)).xyBox(count2,:) = reshape([points([3 4 5],:);addPoint]',1,[]);
%             cuboid.xyBox(i,9:16) = reshape([points(rem(pointIDs(2)+[0 1 2]-1,6)+1,:);addPoint]',1,[]);
%             cuboid.views(i,2) = setdiff(1:3, direct_assign([3 4]));
% %             count3 = rectangle(directIDs(3)).count + 1;
% %             rectangle(directIDs(3)).count = rectangle(directIDs(3)).count + 1;
% %             rectangle(directIDs(3)).xyBox(count3,:) = reshape([points([5 6 1],:);addPoint]',1,[]);
%             cuboid.xyBox(i,17:24) = reshape([points(rem(pointIDs(3)+[0 1 2]-1,6)+1,:);addPoint]',1,[]);
%             cuboid.views(i,3) = setdiff(1:3, direct_assign([5 6]));
%             
% %             cuboid.views(i,:) = directIDs;
% %             cuboid.combo(i,:) = [count1 count2 count3];
%             cuboid.score(i) = score(i);
% %         end
%     end 
% end
% valid = cuboid.score>0.8;
% cuboid.views = cuboid.views(valid,:);
% cuboid.score = cuboid.score(valid);
% cuboid.xyBox = cuboid.xyBox(valid,:);
% % cuboid.combo = cuboid.combo(valid,:);
% % cuboid.xyzBox = zeros(length(cuboid.score),12);
% % for i = 1:length(cuboid.score)
% %     
% % end
% % 


end

