function compObjectHypot( aid )
%COMPOBJECTHYPOT Generate object hypothesis
%   You should have run compRoomHypot first

global config;
bufname = config.bufname;

rot_file_name = [bufname config.IMGLIST(aid).name config.rotateImageSuffix];
rotImg = imread(rot_file_name);

% selective search
SS = getSelectiveSearch(rotImg);
AllCuboid = regionBasedHypothesisFromSS(SS, 1024, 2048);
rectangle = rectDetectionFlexFromSS(SS, 1024, 2048);
AllCuboid_NEW = pruneCuboid(rotImg, AllCuboid);
rectangle_NEW = pruneRectangle(rotImg, rectangle);
% rectangle hypothesis
[ rectangle, rectCuboid2, rectCuboid3 ] = rectangleBasedHypothesis( rotImg );   
% segmentation hypothesis
[ regionCuboid ] = regionBasedHypothesis( rotImg );
% merge up everything
regionCuboid.views = [regionCuboid.views; AllCuboid_NEW.views];
regionCuboid.xyzBox = [regionCuboid.xyzBox; AllCuboid_NEW.xyzBox];
regionCuboid.score = [regionCuboid.score; AllCuboid_NEW.score];
regionCuboid.count = regionCuboid.count + AllCuboid_NEW.count;

for i = 1:6
    rectangle(i).highPixelBox = [rectangle(i).highPixelBox; rectangle_NEW(i).highPixelBox];
    rectangle(i).xyzBox = [rectangle(i).xyzBox; rectangle_NEW(i).xyzBox];
    rectangle(i).score = [rectangle(i).score; rectangle_NEW(i).score];
end

if config.UPDATE_TO_DISK
    parsave([bufname config.objHypoFileUpdate], ...
        'rectangle', rectangle, 'rectCuboid2', rectCuboid2, 'rectCuboid3', rectCuboid3, ...
        'regionCuboid', regionCuboid);
end

[ rect_score ] = getSegScrOfRect( rotImg, rectangle );
CAN_POOL = generateObjHyps( aid, 20000, rect_score);
if config.UPDATE_TO_DISK
    save([bufname config.rectangleScore], 'rect_score');
    save([bufname config.objHypoFile3D],'CAN_POOL');
end
fprintf('finish.\n');
end

