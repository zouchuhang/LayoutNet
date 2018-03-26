function compRoomHypot( aid )
%compRoomHypot compute line segment detection, vanishing point, room layout
%hypothesis, feature like OM, GC, CN are also computed because they rely on
%line segments.
%   To be consistent with previous buffer data, take care of GC

global config
bufname = config.bufname;

%%
Img = imread([config.folderName config.IMGLIST(aid).name '/' config.IMGLIST(aid).name '.jpg']);
Img = im2double(Img);
imgSize = 320;
qError = 0.7;
[ olines, vp, views, edges, panoEdge, score, angle] ...
    = panoEdgeDetection( Img, imgSize, qError); 

Img_small = imresize(Img, [1024 2048]);
[rotImg, R] = rotatePanorama(Img_small, vp(3:-1:1,:));

if config.UPDATE_TO_DISK
    sml_file_name = [bufname config.IMGLIST(aid).name config.smallImageSuffix];
    imwrite(Img_small, sml_file_name);

    rot_file_name = [bufname config.IMGLIST(aid).name config.rotateImageSuffix];
    imwrite(rotImg, rot_file_name);

    parsave([bufname config.vpEstimationFile], 'vp', vp, 'R', R);
end

%%
[ ~, panoOmap ] = computePanoOmap( views, edges, vp );
panoOmap_rot = rotatePanorama(panoOmap, [], R);   
[ ~, wallPanoNormal] = compSurfaceLabel( rotImg );
wallPanoNormal_rot = rotatePanorama(wallPanoNormal, [], R);

orientation_ori = panoOmap; %C.orientation_ori;
surfacelabel = rotatePanorama(wallPanoNormal, inv(R)); %rotatePanorama(C.surfacelabel_ori, inv(R));
hyps = generateHypsB(olines, vp, 3, orientation_ori, surfacelabel);
hyps_rot = rotateHyps( hyps, R);

rotImg = min(max(rotImg,0),1);
colorName_rot = im2cM(double(rotImg));
colorName = rotatePanorama(colorName_rot, inv(R));

if config.UPDATE_TO_DISK
    parsave([bufname config.roomModelFile], ...
        'orientation_ori', panoOmap, 'surfacelabel_ori', wallPanoNormal, ...
        'orientation', panoOmap_rot, 'surfacelabel', wallPanoNormal_rot, ...
        'colorName', colorName, 'colorName_rot', colorName_rot);
    parsave([bufname config.roomHypoFile], ...
        'hypothesis_ori', hyps, 'hypothesis', hyps_rot);
end

end

