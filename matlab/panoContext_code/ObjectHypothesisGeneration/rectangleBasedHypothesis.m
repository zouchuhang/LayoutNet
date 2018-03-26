function [ rectangle, RuleCuboid, RuleCuboidComb ] = rectangleBasedHypothesis( rotImg )
%RECTANGLEBASEDHYPOTHESIS Generate detection based object hypothesis
%   rectangle: rectangle detected in 6 views
%   RuleCuboid: combine two perpendicular rectangles and form cuboids
%       views: id of normal direction of surfaces
%       comb: rect id in corresponding normal id
%   RuleCuboidComb: combine three perpendicular rects to form cuboids

%% wall views
vp = [-1  0  0; ...
       1  0  0; ...
       0 -1  0; ...
       0  1  0; ...
       0  0 -1; ...
       0  0  1];
uv = xyz2uvN(vp);
fov = 160/180*pi;
cutSize = 2000;
rotImg = im2double(rotImg);
[sepSceneSix] = separatePano( rotImg, fov, uv(:,1), uv(:,2), cutSize);

%% load detector
% load './rectangleDetector/finalModel.mat'
% rectModel = allModel{4}; % after 2nd negative mining

%% rectangle detect	
detconfig;
load(config.modelfile);

% try
%     load([bufname '/rectDetResult.mat']);
% catch
    rectangle = repmat(struct('highPixelBox',[]), 6, 1);
    for vid = 1:6
        img = sepSceneSix(vid).img;
        [ highPixelBox, highFeatureBox, feature] = rectDetectionFlex( img.*255, rectModel, config );
        rectangle(vid).highPixelBox = highPixelBox;
        rectangle(vid).fov = sepSceneSix(vid).fov;
        rectangle(vid).vx = sepSceneSix(vid).vx;
        rectangle(vid).vy = sepSceneSix(vid).vy;
        rectangle(vid).cutSize = size(sepSceneSix(vid).img,1);
    end
%     save( [bufname '/rectDetResult.mat'], 'rectangle', 'config');
% end

%% combine rectangle to cuboids
% try
%     load( [bufname '/rectangleBasedCuboid.mat'] );
% catch
    %% convert pixel box to uv coordinate
%     resizeImg = imresize(rotImg, [512 1024]);
    for vid = 1:6
        highPixelBox = rectangle(vid).highPixelBox;
        highPixelBox = highPixelBox(highPixelBox(:,end)>config.newRepThresh,: );

        viewSize = size(sepSceneSix(vid).img, 1);
        R = viewSize/2/tan(sepSceneSix(vid).fov/2);
        pt1 = [highPixelBox(:,1) highPixelBox(:,2)];
        pt2 = [highPixelBox(:,3) highPixelBox(:,2)];
        pt3 = [highPixelBox(:,3) highPixelBox(:,4)];
        pt4 = [highPixelBox(:,1) highPixelBox(:,4)];

        xyz0 = R*uv2xyzN([sepSceneSix(vid).vx sepSceneSix(vid).vy], 1);
        [ xyz1, ~ ] = xyzFromView( pt1, xyz0, viewSize, viewSize );
        [ xyz2, ~ ] = xyzFromView( pt2, xyz0, viewSize, viewSize );
        [ xyz3, ~ ] = xyzFromView( pt3, xyz0, viewSize, viewSize );
        [ xyz4, ~ ] = xyzFromView( pt4, xyz0, viewSize, viewSize );

        rectangle(vid).xyzBox = [xyz1 xyz2 xyz3 xyz4];
        rectangle(vid).score = highPixelBox(:,5);
    end

    %% combine 2 rectangles
    combrules;
    for rid = 1:8
        rule = CubRule(rid);
        viewID1 = rule.suf(1);
        viewID2 = rule.suf(2);
        view1point1 = find(GeoRule(viewID1).pointID==rule.cont(1));
        view1point2 = find(GeoRule(viewID1).pointID==rule.cont(2));
        view2point1 = find(GeoRule(viewID2).pointID==rule.cont(1));
        view2point2 = find(GeoRule(viewID2).pointID==rule.cont(2));

        lhsBox = rectangle(viewID1).xyzBox;
        lhsScore = rectangle(viewID1).score;
        rhsBox = rectangle(viewID2).xyzBox;
        rhsScore = rectangle(viewID2).score;

        lhsPts = lhsBox(:,[view1point1*3-2:view1point1*3 view1point2*3-2:view1point2*3]);
        lhsNormal = cross(lhsPts(:,1:3), lhsPts(:,4:6), 2);
        lhsNormal = lhsNormal ./ repmat( sqrt(sum(lhsNormal.^2,2)), 1, 3);
        lhsAngLen = acos( sum(lhsPts(:,1:3).*lhsPts(:,4:6),2));

        rhsPts = rhsBox(:,[view2point1*3-2:view2point1*3 view2point2*3-2:view2point2*3]);
        rhsNormal = cross(rhsPts(:,1:3), rhsPts(:,4:6), 2);
        rhsNormal = rhsNormal ./ repmat( sqrt(sum(rhsNormal.^2,2)), 1, 3);
        rhsAngLen = acos( sum(rhsPts(:,1:3).*rhsPts(:,4:6),2));

        distMatrix = zeros( size(rhsPts,1), size(lhsPts,1)); 
        normMatrix = zeros( size(rhsPts,1), size(lhsPts,1));
        distMat1 = zeros( size(rhsPts,1), size(lhsPts,1)); 
        distMat2 = zeros( size(rhsPts,1), size(lhsPts,1)); 
        for i = 1:size(lhsPts,1)
            angle = acos(abs(sum(( repmat(lhsNormal(i,1:3), size(rhsPts,1), 1).*rhsNormal(:,1:3)), 2)));        
            dist1 = acos(sum((repmat(lhsPts(i,1:3), size(rhsPts,1), 1).*rhsPts(:,1:3)), 2));
            dist2 = acos(sum((repmat(lhsPts(i,4:6), size(rhsPts,1), 1).*rhsPts(:,4:6)), 2));

            distMat1(:,i) = dist1;
            distMat2(:,i) = dist2;
            normMatrix(:,i) = angle;
            distMatrix(:,i) = min(rhsAngLen,repmat(lhsAngLen(i), size(rhsPts,1), 1));
    %         distMatrix(:,i) = (max(dist1,dist2))./(min(rhsAngLen,repmat(lhsAngLen(i), size(rhsPts,1), 1)));

        end

        normValid = normMatrix<0.04;
        distValid = max(distMat1, distMat2)./min(distMatrix, 1)<0.15;

        [matchRhsID, matchLhsID] = find(normValid & distValid);
    %     RuleCuboid(rid).matchLhsID = matchLhsID;
    %     RuleCuboid(rid).matchRhsID = matchRhsID;
        RuleCuboid(rid).combo = [matchLhsID matchRhsID zeros(length(matchLhsID),1)];
        RuleCuboid(rid).views = [CubRule(rid).suf 0];
        RuleCuboid(rid).score = lhsScore(matchLhsID) + rhsScore(matchRhsID);
        RuleCuboid(rid).lhsAngLen = lhsAngLen(matchLhsID);
        RuleCuboid(rid).rhsAngLen = rhsAngLen(matchRhsID);

        I = sub2ind([size(rhsPts,1) size(lhsPts,1)], matchRhsID, matchLhsID);

        RuleCuboid(rid).dist0 = normMatrix(I);
        RuleCuboid(rid).dist1 = max(distMat1(I), distMat2(I));
        RuleCuboid(rid).dist2 = min(distMatrix(I), 1);
    end

    %% combine 3 rectangles
    RULESURF = vertcat(CubRule.suf);

    for rid = 1:4
        views = sort([CubRule(rid).suf CubRule(rid).check], 'ascend');
        rules(1) = find(RULESURF(:,1)==views(1) & RULESURF(:,2)==views(2) | (RULESURF(:,1)==views(2) & RULESURF(:,2)==views(1)));
        rules(2) = find(RULESURF(:,1)==views(2) & RULESURF(:,2)==views(3) | (RULESURF(:,1)==views(3) & RULESURF(:,2)==views(2)));
        rules(3) = find(RULESURF(:,1)==views(3) & RULESURF(:,2)==views(1) | (RULESURF(:,1)==views(1) & RULESURF(:,2)==views(3)));

        RuleCuboidComb(rid).views = views;
        RuleCuboidComb(rid).rules = rules;

        for i = 1:3
            if views(i) == CubRule(rules(i)).suf(1)
                rectList(i).list1 = RuleCuboid(rules(i)).combo(:,1);
                rectList(i).list2 = RuleCuboid(rules(i)).combo(:,2);
            else
                rectList(i).list1 = RuleCuboid(rules(i)).combo(:,2);
                rectList(i).list2 = RuleCuboid(rules(i)).combo(:,1);
            end
        end
        % from 1st to 2nd   
        combo = zeros(1000,3);
        count = 0;

        for cid1 = 1:length(rectList(1).list1)
            validID1 = rectList(1).list1(cid1);
            validID2 = rectList(1).list2(cid1);
            validID3 = rectList(2).list2( rectList(2).list1 == validID2);
            validID3_plus = rectList(3).list1(rectList(3).list2 == validID1);
            C = intersect(validID3, validID3_plus);
            for i = 1:length(C)
                count = count + 1;
                combo(count,:) = [validID1 validID2 C(i)];
            end
        end
        combo = combo(1:count,:);
        RuleCuboidComb(rid).combo = combo;
    end

%     save( [bufname '/rectangleBasedCuboid.mat'], 'rectangle', 'RuleCuboid', 'RuleCuboidComb');
% end

end

