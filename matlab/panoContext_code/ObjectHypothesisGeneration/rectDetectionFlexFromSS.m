function [ rectangle ] = rectDetectionFlexFromSS( SS1, imgH, imgW )
%RECTANGLEDETECTIONFROMSS Convert selective search to rectangles
%   SS1: selective search hierarchy grouping of segmentation
%   rectangle:
%       xyzBox: write 4 corners to a one dimensional vector
%       score: intersection/union with segmentation


% colorTypes = {'Hsv', 'Lab', 'RGI', 'H', 'Intensity'};
vp = [-1  0  0; ...
       1  0  0; ...
       0 -1  0; ...
       0  1  0; ...
       0  0 -1; ...
       0  0  1];
uv = xyz2uvN(vp);
fov = 160/180*pi;
cutSize = 500;
% [sepSceneSix] = separatePano( rotImg, fov, uv(:,1), uv(:,2), cutSize);

rectangle = repmat(struct('highPixelBox',zeros(2000,5),'count',0), 6, 1);
min_size = 2500;
max_size = imgH*imgW/2;
big_map = zeros(imgH, imgW);
[gridX, gridY] = meshgrid(1:cutSize, 1:cutSize);
SE = strel('disk', 5, 4);
SS = strel('disk', 2, 4);


for colorid = 1:length(SS1)
    fprintf('colorID: %d/%d\n', colorid, length(SS1));
    hBlobs1 = RecreateBlobHierarchyIndIm( SS1(colorid).blobIndIm, SS1(colorid).blobBoxes, SS1(colorid).hierarchy{1});
    hBlobs2 = RecreateBlobHierarchyIndIm( SS1(colorid).blobIndIm, SS1(colorid).blobBoxes, SS1(colorid).hierarchy{2});
    hBlobs = [hBlobs1; hBlobs2];
    
%     hBlobs_valid = false(length(hBlobs),1);
%     for hid = 1:length(hBlobs)
%         if hBlobs{hid}.size>min_size && hBlobs{hid}.size<max_size
%             hBlobs_valid(hid) = true;
%         end
%     end
    [ hBlobs_valid, hBlobs_centroid ] = ...
        regionShapeValidation( hBlobs, 1024, 2048, 2500, 3, 0.65, false);
       
    for bid = find(hBlobs_valid')
        fprintf('BlobID: %d/%d\n', bid, length(hBlobs_valid));
        rect = hBlobs{bid}.rect;
        map = big_map;
        map(rect(1):rect(3),rect(2):rect(4)) = hBlobs{bid}.mask;
        rect_ctr_coord = [(rect(2)+rect(4))/2 (rect(1)+rect(3))/2];
        rect_ctr_xyz = uv2xyzN(coords2uv(rect_ctr_coord, 2048, 1024));
        test_vp_ids = find( sum(vp.*repmat(rect_ctr_xyz, 6, 1), 2)>0);
        
        [sepScene] = separatePano( map, fov, uv(test_vp_ids,1), uv(test_vp_ids,2), cutSize);
        for sid = 1:length(sepScene)
%             binaryMap = imdilate( sepScene(sid).img, SE);
%             binaryMap = imerode( binaryMap, SE);
%             binaryMap = imerode( binaryMap, SS);
%             binaryMap = imdilate( binaryMap, SS);
%             binaryMap = imfill(binaryMap, 'hole');
            binaryMap = sepScene(sid).img;
            indX = gridX(binaryMap(:)>0);
            indY = gridY(binaryMap(:)>0);
%             indX = gridX(sepScene(sid).img(:)>0);
%             indY = gridY(sepScene(sid).img(:)>0);
            if ~isempty(indX) && ~isempty(indY) && min(indX)~=1 && min(indY)~=1 && max(indX)~=cutSize && max(indY)~=cutSize
                loc = rectangle(test_vp_ids(sid)).count + 1;
                scr = sum(sepScene(sid).img(:)>0)/(max(indY)-min(indY))/(max(indX)-min(indX));
                rectangle(test_vp_ids(sid)).highPixelBox(loc,:) = ...
                    [min(indX) min(indY) max(indX) max(indY) scr];
                rectangle(test_vp_ids(sid)).count = loc;
            end
        end
    end   
end

for vid = 1:6
    rectangle(vid).highPixelBox = rectangle(vid).highPixelBox(1:rectangle(vid).count,:);
    rectangle(vid).fov = fov;
    rectangle(vid).vx = uv(vid,1);
    rectangle(vid).vy = uv(vid,2);
    rectangle(vid).cutSize = cutSize;
end

for vid = 1:6
    highPixelBox = rectangle(vid).highPixelBox;

    R = cutSize/2/tan(fov/2);
    pt1 = [highPixelBox(:,1) highPixelBox(:,2)];
    pt2 = [highPixelBox(:,3) highPixelBox(:,2)];
    pt3 = [highPixelBox(:,3) highPixelBox(:,4)];
    pt4 = [highPixelBox(:,1) highPixelBox(:,4)];

    xyz0 = R*uv2xyzN(uv(vid,:), 1);
    [ xyz1, ~ ] = xyzFromView( pt1, xyz0, cutSize, cutSize );
    [ xyz2, ~ ] = xyzFromView( pt2, xyz0, cutSize, cutSize );
    [ xyz3, ~ ] = xyzFromView( pt3, xyz0, cutSize, cutSize );
    [ xyz4, ~ ] = xyzFromView( pt4, xyz0, cutSize, cutSize );

    rectangle(vid).xyzBox = [xyz1 xyz2 xyz3 xyz4];
    rectangle(vid).score = highPixelBox(:,5);
end

end

