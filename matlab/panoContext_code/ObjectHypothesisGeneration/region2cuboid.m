function [ maxHypot ] = region2cuboid( binaryMap, vanishing )
%REGION2CUBOID Summary of this function goes here
%   Detailed explanation goes here

binaryMap = imresize(binaryMap, 0.5);
binaryMap = binaryMap>0.5;
vanishing = vanishing .* 0.5;

%% process the mask: dilate and erode to smooth and kill seam, then fill hole
% SE = strel('disk', 5, 4);
% binaryMap = imdilate( binaryMap, SE);
% binaryMap = imerode( binaryMap, SE);
% binaryMap = imfill(binaryMap, 'hole');
[imgH, imgW, ~] = size(binaryMap);

% get region of interesting to accelerate all computation
XSum = find(sum(binaryMap, 1)>0);
YSum = find(sum(binaryMap, 2)>0);
tlX = max(min(XSum)-10,1);    boxW = min(max(XSum)+10,imgW)-tlX+1;
tlY = max(min(YSum)-10,1);    boxH = min(max(YSum)+10,imgH)-tlY+1;

binaryMap = binaryMap(tlY:tlY+boxH-1, tlX:tlX+boxW-1);
vanishing(:,1) = vanishing(:,1) - tlX + 1;
vanishing(:,2) = vanishing(:,2) - tlY + 1;

%% fprintf('compute angle for all pixels\n');
% tic;
[imgH, imgW, ~] = size(binaryMap);
[gridX, gridY] = meshgrid( 1:imgW, 1:imgH);
angleMap = zeros(imgH, imgW, 3);
for i = 1:3
    R = atan((gridY-vanishing(i,2))./(gridX-vanishing(i,1)+0.00001));
    R(R<0 & gridY-vanishing(i,2)>0) = R(R<0 & gridY-vanishing(i,2)>0)+pi;
    R(R>0 & gridY-vanishing(i,2)<0) = R(R>0 & gridY-vanishing(i,2)<0)-pi;
    angleMap(:,:,i) = R;
end
% toc;

%% fprintf('Define boundary and interior region\n');
boundary = edge(binaryMap, 'canny', [0.1 0.2], 1);
distMap = bwdist(boundary);
BoundaryRegion = distMap<5 & binaryMap;
maxDistance = max(distMap(binaryMap));
InsideRegion = binaryMap & distMap>maxDistance/2;

%% fprintf('define possible test range\n')
angleTest = repmat(struct('angle',[],'lines',[],'vanish',[],'direct',[]), 6, 1);
for i = 1:3   
    map = angleMap(:,:,i);
    BR = map(BoundaryRegion(:));
    minBR = min(BR(:));
    maxBR = max(BR(:));
    IR = map(InsideRegion(:));
    minIR = min(IR(:));
    maxIR = max(IR(:));
    
    if maxBR-minBR<pi;
        rangeBR = [minBR maxBR];
    else
        rangeBR = [min(BR(BR>0)) max(BR(BR<0))];
    end
    if maxIR-minIR<pi;
        rangeIR = [minIR maxIR];
    else
        rangeIR = [min(IR(IR>0)) max(IR(IR<0))];
    end
    % LHS
    angleTest(2*i-1).angle = uniformSample( rangeBR(1), rangeIR(1), 20);    
    angleTest(2*i-1).lines = [repmat(vanishing(i,:), 20, 1)-1000*[cos(angleTest(2*i-1).angle)' sin(angleTest(2*i-1).angle)'] ...
        repmat(vanishing(i,:), 20, 1)+1000*[cos(angleTest(2*i-1).angle)' sin(angleTest(2*i-1).angle)']];
    angleTest(2*i-1).vanish = vanishing(i,:);
    angleTest(2*i-1).direct = [cos(angleTest(2*i-1).angle)' sin(angleTest(2*i-1).angle)' zeros(20,1)];
    %RHS
    angleTest(2*i).angle = uniformSample( rangeIR(2), rangeBR(2), 20);   
    angleTest(2*i).lines = [repmat(vanishing(i,:), 20, 1)-1000*[cos(angleTest(2*i).angle)' sin(angleTest(2*i).angle)'] ...
    repmat(vanishing(i,:), 20, 1)+1000*[cos(angleTest(2*i).angle)' sin(angleTest(2*i).angle)']];
    angleTest(2*i).vanish = vanishing(i,:);
    angleTest(2*i).direct = [cos(angleTest(2*i).angle)' sin(angleTest(2*i).angle)' zeros(20,1)];
end
% figure; clf;
% imshow(binaryMap); hold on;
% color = [1 0 0; 0.5 0 0; 0 1 0; 0 0.5 0; 0 0 1; 0 0 0.5];
% for i = 1:6
%     for j = 1:length(angleTest(i).angle)
%         L = [angleTest(i).vanish-1000*[cos(angleTest(i).angle(j)) sin(angleTest(i).angle(j))] ...
%              angleTest(i).vanish+1000*[cos(angleTest(i).angle(j)) sin(angleTest(i).angle(j))]];
%         line(L([1 3]), L([2 4]), 'Color', color(i,:), 'LineWidth', 2);
%     end
% end

%% fprintf('new ransac: two stages');
% first stage
comb = [
     1     3; ...
     1     4; ...
     1     5; ...
     1     6; ...
     2     3; ...
     2     4; ...
     2     5; ...
     2     6; ...
     3     5; ...
     3     6; ...
     4     5; ...
     4     6];
interCheck = [3 3 2 2 3 3 2 2 1 1 1 1];
vanishCheck1 = vertcat(angleTest(interCheck.*2-1).vanish);
vanishCheck2 = vertcat(angleTest(interCheck.*2).vanish);
maxTrial = 10000;
maxValidNum = 250;
Ls = zeros(6,4);
DT = zeros(6,3);
valid_points = repmat(struct('direct',zeros(3,1),'points',zeros(6,2), 'type', 0, 'startID', zeros(1,3), 'addPoint', zeros(1,2)),maxValidNum,1);
num = 0;
iter = 0;
while iter<maxTrial && num<maxValidNum
    iter = iter + 1;
    for i = 1:6
        I = randi(20,1);
        Ls(i,:) = angleTest(i).lines(I,:);
        DT(i,:) = angleTest(i).direct(I,:);
    end
    intersect = intersectPoint( Ls(comb(:,1),:), Ls(comb(:,2),:));
    cross1 = cross([intersect-vanishCheck1 zeros(12,1)], DT(interCheck.*2-1,:), 2);
    cross2 = cross([intersect-vanishCheck2 zeros(12,1)], DT(interCheck.*2,:), 2);
    
    valid = cross1(:,3)<0 & cross2(:,3)>0; % take the 6 points
    if sum(valid)~=6
        continue;
    else
        intersect = intersect(valid,:); % check convex
        K = convhull(intersect);
        if size(K,1)~=7
            continue;
        else
            points = intersect(K(1:6),:);
            %[x2,y2] = poly2cw(points(:,1), points(:,2)); 
		[x2,y2] = poly2cw_N(points(:,1), points(:,2));
            points = [x2 y2];
        end       
    end
    
    lines = points([2 3 4 5 6 1],:) - points([1 2 3 4 5 6],:);
    direct = lines./repmat(sqrt(sum(lines.^2,2)), 1, 2);
    midpoints = (points([2 3 4 5 6 1],:) + points([1 2 3 4 5 6],:))/2;
    A = zeros(6,3);
    for j = 1:3
        vanline = [midpoints(:,1)-vanishing(j,1) midpoints(:,2)-vanishing(j,2)];
        vandire = vanline./repmat(sqrt(sum(vanline.^2,2)), 1, 2);
        A(:,j) = dot(direct, vandire, 2);
    end
    [~,direct_assign] = max( abs(A), [], 2);
    type_indicate = direct_assign(1:3)==direct_assign(4:6);
    
    
    if sum(type_indicate)==3 && all(ismember([1 2 3], direct_assign)) % type 3
        I = find(direct_assign==3);
%         if numel(I)==0
%             fprintf('');
%         end
        verticalID = rem([I(1)+2 I(1)+5]-1,6) + 1;
        if norm(points(verticalID(1),:)-vanishing(3,:))<norm(points(verticalID(2),:)-vanishing(3,:))
            pointIDs = rem([verticalID(1) verticalID(1)+2 verticalID(1)+4]-1,6)+1;
        else
            pointIDs = rem([verticalID(2) verticalID(2)+2 verticalID(2)+4]-1,6)+1;
        end
        directIDs = direct_assign(rem(pointIDs-1+1,6)+1);
        lines = [points(pointIDs(1),:) vanishing(directIDs(1),:); ...
                 points(pointIDs(2),:) vanishing(directIDs(2),:); ...
                 points(pointIDs(3),:) vanishing(directIDs(3),:)];
        intersect = intersectPoint(lines([1 2 3],:), lines([2 3 1],:));
        if sum(sum((intersect([1 2 3],:)-intersect([2 3 1],:)).^2,2))<27 % pass the intersection test, 3 pixel error
            num = num + 1;
            valid_points(num).points = points;
            valid_points(num).type = 3;
            valid_points(num).startID = pointIDs;
            valid_points(num).addPoint = sum(intersect,1)./3;
            direct1 = setdiff(1:3,direct_assign(rem(([pointIDs(1):pointIDs(1)+1]-1),6)+1));
            direct2 = setdiff(1:3,direct_assign(rem(([pointIDs(2):pointIDs(2)+1]-1),6)+1));
            direct3 = setdiff(1:3,direct_assign(rem(([pointIDs(3):pointIDs(3)+1]-1),6)+1));
            valid_points(num).direct = [direct1 direct2 direct3];
        end
        
    elseif sum(type_indicate)==1  && all(ismember([1 2 3], direct_assign)) % type 1
        I = find(type_indicate);
        if ~(direct_assign(rem(I-1+1,6)+1)==direct_assign(rem(I-1+5,6)+1) ...
                && direct_assign(rem(I-1+2,6)+1)==direct_assign(rem(I-1+4,6)+1))
            continue;
        end
        
        pointID1 = rem((I-1)+2,6)+1;
        pointID2 = rem((pointID1-1)+3,6)+1;
        testdir = points(pointID1,:)-points(pointID2,:);
        testdir = testdir ./ norm(testdir);
        refedir = vanishing(direct_assign(I),:) - (points(pointID1,:)+points(pointID2,:))/2;
        refenorm = norm(refedir);
        
        distance = refenorm * sqrt(1-abs(dot(testdir, refedir./refenorm, 2))^2);
        if distance<16 % pass the parallel test, ~8 degree error
            num = num + 1;
            valid_points(num).points = points;
            valid_points(num).type = 1;
            valid_points(num).startID(1:2) = [pointID1 pointID2];
            direct1 = setdiff(1:3,direct_assign(rem(([pointID1:pointID1+2]-1),6)+1));
            direct2 = setdiff(1:3,direct_assign(rem(([pointID2:pointID2+2]-1),6)+1));
            valid_points(num).direct = [direct1 direct2 0];            
        end
    end

end

% check consistency

if num==0
    maxHypot = struct('direct',zeros(3,1),'points',zeros(6,2), 'type', 0, 'startID', zeros(1,3), 'addPoint', zeros(1,2));
    maxHypot.score = 0;
else
    scores = zeros(num,1);
    for tid = 1:num
        BW = poly2mask( valid_points(tid).points(:,1),  valid_points(tid).points(:,2), imgH, imgW);
        IS = BW & binaryMap; UN = BW | binaryMap;
        SC = sum(IS(:))/sum(UN(:));
        scores(tid) = SC;
    end
    [B,I] = max(scores);
    maxHypot = valid_points(I);
    points = maxHypot.points;
    addPoint = maxHypot.addPoint;
    maxHypot.points = [points(:,1)+tlX-1 points(:,2)+tlY-1] .* 2;
    maxHypot.addPoint = [addPoint(1)+tlX-1 addPoint(2)+tlY-1] .* 2;
    maxHypot.score = B;
end
%% fprintf('old ransac\n');
% % connectRule = [1 6 3 2 5 4 1];
% comb = [
%      1     3; ...
%      1     4; ...
%      1     5; ...
%      1     6; ...
%      2     3; ...
%      2     4; ...
%      2     5; ...
%      2     6; ...
%      3     5; ...
%      3     6; ...
%      4     5; ...
%      4     6];
% interCheck = [3 3 2 2 3 3 2 2 1 1 1 1];
% vanishCheck1 = vertcat(angleTest(interCheck.*2-1).vanish);
% vanishCheck2 = vertcat(angleTest(interCheck.*2).vanish);
% maxIter = 500;
% maxScore = 0;
% iter  = 0;
% Ls = zeros(6,4);
% DT = zeros(6,3);
% while iter<maxIter && maxScore<0.85
%     for i = 1:6
%         I = randi(20,1);
%         Ls(i,:) = angleTest(i).lines(I,:);
%         DT(i,:) = angleTest(i).direct(I,:);
%     end
%     intersect = intersectPoint( Ls(comb(:,1),:), Ls(comb(:,2),:));
% 
% %     figure(100);  imshow(binaryMap); hold on
% %     for i = 1:6
% %         line(Ls(i,[1 3]), Ls(i,[2 4]), 'Color', color(i,:), 'LineWidth', 2);
% %     end 
% 
%     cross1 = cross([intersect-vanishCheck1 zeros(12,1)], DT(interCheck.*2-1,:), 2);
%     cross2 = cross([intersect-vanishCheck2 zeros(12,1)], DT(interCheck.*2,:), 2);
%     valid = cross1(:,3)<0 & cross2(:,3)>0;
%     if sum(valid)~=6
%         continue;
%     else
%         intersect = intersect(valid,:);
%         K = convhull(intersect);
%         if size(K,1)~=7
%             continue;
%         else
%             intersect = intersect(K(1:6),:);
%         end       
%     end
%     
% %     IDs = zeros(12,1);
% %     for i = 1:6
% %         validID = find(comb(:,1)==i | comb(:,2)==i);
% %         points = intersect( comb(:,1)==i | comb(:,2)==i,:);
% %         dist = sum((points- repmat(angleTest(i).vanish, 4, 1)).^2,2);
% %         [~, I] = sort(dist, 'ascend');
% %         IDs(i) = validID(I(2));
% %         IDs(i+6) = validID(I(3));
% %     end
% %     I = unique(IDs);
% %     if length(I)==6
% %         K = convhull(intersect(I,:));
% %         if size(K,1)~=7
% %             continue;
% %         else
% %             intersect = intersect(I(K(1:6)),:);
% %         end
% %     else
% %         continue;
% %     end
%   
% %     
% %     for i =1:6
% %         plot(intersect(i,1), intersect(i,2), '--rs','LineWidth',2,...
% %                 'MarkerEdgeColor','k',...
% %                 'MarkerFaceColor','b',...
% %                 'MarkerSize',10);
% %     end
%     
%     BW = poly2mask( intersect(:,1), intersect(:,2), imgH, imgW);
% %     BW = poly2mask( intersect(:,1)-tlX+1, intersect(:,2)-tlY+1, boxH, boxW);
%     IS = BW & binaryMap; UN = BW | binaryMap;
%     SC = sum(IS(:))/sum(UN(:));
%     
%     if SC>maxScore
%         maxScore = SC;
%         maxHypot = intersect;
% %         maxSource = Ls;
%     end 
%     iter = iter + 1;
% end
% % maxHypot = [maxHypot(:,1) maxHypot(:,2)] .* 2;
% maxHypot = [maxHypot(:,1)+tlX-1 maxHypot(:,2)+tlY-1] .* 2;
% % maxHypot = [maxHypot(:,1)+tlX maxHypot(:,2)+tlY]

end

function X = uniformSample(lowBound, highBound, num)
if highBound<lowBound
    highBound = highBound + 2*pi;
end
X = linspace(lowBound, highBound, num);
X = normValue( X, -pi, pi);
end

function [nv] = normValue( v, minBound, maxBound)
nv = rem((v-minBound)+2*(maxBound-minBound), maxBound-minBound) + minBound;
end
