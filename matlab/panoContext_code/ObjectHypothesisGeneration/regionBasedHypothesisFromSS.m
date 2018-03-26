function [ AllCuboid ] = regionBasedHypothesisFromSS( SS1, imgH, imgW )
%REGIONBASEDHYPOTHESISFROMSS Convert selective search result to cuboid
%   SS1: selective search result: hierarchy grouping of over-segment
%   AllCuboid:
%       xyzbox: write coordinates of 4 corners to a vector
%       views: visible normal direction of surfaces
%       score: intersection/union with segmentation

% [imgH, imgW, ~] = size(rotImg);
vanishingPoint = [-1 0 0; 1 0 0; 0 -1 0; 0 1 0;0 0 -1; 0 0 1];
viewID = [1 2 3 4 5 6];
colorTypes = {'Hsv', 'Lab', 'RGI', 'H', 'Intensity'};

%% extract ss cuboid candidate
AllCuboid = struct('views',zeros(2000,3),'xyzBox',zeros(2000,36),'score',zeros(2000,1),'count',0);
minAngleDist = pi/18;

for colorid = 1:length(colorTypes)
    hBlobs1 = RecreateBlobHierarchyIndIm( SS1(colorid).blobIndIm, SS1(colorid).blobBoxes, SS1(colorid).hierarchy{1});
    hBlobs2 = RecreateBlobHierarchyIndIm( SS1(colorid).blobIndIm, SS1(colorid).blobBoxes, SS1(colorid).hierarchy{2});
    hBlobs = [hBlobs1; hBlobs2];
    [ hBlobs_valid, hBlobs_centroid ] = regionShapeValidation( hBlobs, 1024, 2048, 2500, 3, 0.65, true );
    
    hBlobs_ViewPoint = coords2uv( hBlobs_centroid, imgW, imgH );
    
    yValid = hBlobs_ViewPoint(:,2)/(minAngleDist);
    yValid(abs(yValid)>=1) = 0;
    hBlobs_ViewPoint(yValid~=0,2) = minAngleDist*sign(yValid(yValid~=0));
    
    for j = 1:size(hBlobs_ViewPoint,1)
        xClose = hBlobs_ViewPoint(j,1) - [-pi -pi/2 0 pi/2 pi];
        [B,I] = min(abs(xClose));
        if B<minAngleDist
            hBlobs_ViewPoint(j,1) = -pi + (I-1)*pi/2 + sign(xClose(I))*minAngleDist;
        end
    end
    
    fov = 2/3*pi;
    fovSize = 640;
    big_map = zeros(imgH, imgW);
    
    for bid = find(hBlobs_valid')
        rect = hBlobs{bid}.rect;
        map = big_map;
        map(rect(1):rect(3),rect(2):rect(4)) = hBlobs{bid}.mask;
        warped_image = imgLookAt(map, hBlobs_ViewPoint(bid,1), hBlobs_ViewPoint(bid,2), fovSize, fov );
        
        fprintf('%d/%d, ', bid, length(hBlobs));
        vp = uv2xyzN(hBlobs_ViewPoint(bid,:));
	fprintf('1');
        [out2D, valid, division] = projectPoint2SeparateView( vanishingPoint, vp, fov, fovSize );
	fprintf('2');
        vanishing = out2D(valid,:);
	fprintf('3');
        subViewID = viewID(valid);
	fprintf('4');
        [ cuboid ] = getSegHypothesisID( warped_image, 200, vanishing, 1);        

	fprintf('Get cuboid\n');
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
    
end

AllCuboid.views = AllCuboid.views(1:AllCuboid.count,:);
AllCuboid.xyzBox = AllCuboid.xyzBox(1:AllCuboid.count,:);
AllCuboid.score = AllCuboid.score(1:AllCuboid.count);

end

