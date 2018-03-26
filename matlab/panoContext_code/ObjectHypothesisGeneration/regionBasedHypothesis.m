function [ AllCuboid ] = regionBasedHypothesis( rotImg )
%REGIONBASEDHYPOTHESIS Generate segmentation based hypothesis by using different
%segmentation parameters
%   Detailed explanation goes here
%%
vanishingPoint = [-1 0 0; 1 0 0; 0 -1 0; 0 1 0;0 0 -1; 0 0 1];
viewID = [1 2 3 4 5 6];

AllCuboid = struct('views',zeros(2000,3),'xyzBox',zeros(2000,36),'score',zeros(2000,1),'count',0);

thresholds = [200 500 800 1200 2000 3000];
minAngleDist = pi/18;
for i = 1:6
    AllViewPt(i).labelMap = gbPanoSegment( im2uint8(rotImg), 0.5, thresholds(i), 50);
%     figure; imshow(AllViewPt(i).labelMap,[]);
    AllViewPt(i).thresh = thresholds(i);
    [viewpoint, labelID, segCentroid] = getValidViewPoint(AllViewPt(i).labelMap, 2500);
    
    yValid = viewpoint(:,2)/(minAngleDist);
    yValid(abs(yValid)>=1) = 0;
    viewpoint(yValid~=0,2) = minAngleDist*sign(yValid(yValid~=0));
    
    for j = 1:size(viewpoint,1)
        xClose = viewpoint(j,1) - [-pi -pi/2 0 pi/2 pi];
        [B,I] = min(abs(xClose));
        if B<minAngleDist
            viewpoint(j,1) = -pi + (I-1)*pi/2 + sign(xClose(I))*minAngleDist;
        end
    end
    
    AllViewPt(i).viewpoint = viewpoint;
    AllViewPt(i).segCentroid = segCentroid;
    AllViewPt(i).labelID = labelID;
    
    fov = 2/3*pi;
    fovSize = 640;
    if size(viewpoint,1) == 0
        continue;
    end
    sepSceneSeg = separatePano(AllViewPt(i).labelMap, fov, viewpoint(:,1), viewpoint(:,2), fovSize);
    
    for j = 1:length(sepSceneSeg)
        fprintf('%d\n', labelID(j));
        vp = uv2xyzN([sepSceneSeg(j).vx sepSceneSeg(j).vy]);
        [out2D, valid, division] = projectPoint2SeparateView( vanishingPoint, vp, fov, fovSize );
        vanishing = out2D(valid,:);
        subViewID = viewID(valid);
        [ cuboid ] = getSegHypothesisID( sepSceneSeg(j).img, 200, vanishing, labelID(j) );
        
        for k = 1:length(cuboid)
            ID = AllCuboid.count + k;
            AllCuboid.score(ID) = cuboid(k).score;
            if cuboid(k).type==1
                points = cuboid(k).points;
                startID = cuboid(k).startID;
                surf1 = points(rem(([startID(1):startID(1)+3]-1),6)+1,:);
                surf2 = points(rem(([startID(2):startID(2)+3]-1),6)+1,:);
                [ out3DNorm1, ~ ] = projectPointFromSeparateView(  surf1, vp, fov, fovSize );
                [ out3DNorm2, ~ ] = projectPointFromSeparateView(  surf2, vp, fov, fovSize );               
                AllCuboid.xyzBox(ID,1:24) = reshape([out3DNorm1;out3DNorm2]', 1, []);
                AllCuboid.views(ID,1:2) = subViewID(cuboid(k).direct(1:2));
            elseif cuboid(k).type==3
                points = cuboid(k).points;
                startID = cuboid(k).startID;
                surf1 = [points(rem(([startID(1):startID(1)+2]-1),6)+1,:);cuboid(k).addPoint];
                surf2 = [points(rem(([startID(2):startID(2)+2]-1),6)+1,:);cuboid(k).addPoint];
                surf3 = [points(rem(([startID(3):startID(3)+2]-1),6)+1,:);cuboid(k).addPoint];
                [ out3DNorm1, ~ ] = projectPointFromSeparateView(  surf1, vp, fov, fovSize );
                [ out3DNorm2, ~ ] = projectPointFromSeparateView(  surf2, vp, fov, fovSize );
                [ out3DNorm3, ~ ] = projectPointFromSeparateView(  surf3, vp, fov, fovSize );
                AllCuboid.xyzBox(ID,1:36) = reshape([out3DNorm1;out3DNorm2;out3DNorm3]', 1, []);
                AllCuboid.views(ID,:) = subViewID(cuboid(k).direct);
            end
                       
        end
        AllCuboid.count = AllCuboid.count + length(cuboid);
    end
%     figure; viewRegionHypothesis(rotImg, AllCuboid);
end

AllCuboid.views = AllCuboid.views(1:AllCuboid.count,:);
AllCuboid.xyzBox = AllCuboid.xyzBox(1:AllCuboid.count,:);
AllCuboid.score = AllCuboid.score(1:AllCuboid.count);

end

