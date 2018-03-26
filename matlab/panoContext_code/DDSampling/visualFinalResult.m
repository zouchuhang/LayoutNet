global config;
load(fullfile(config.bufname, config.objHypoFile3D));
%%
load(fullfile(config.bufname, config.globalResultFile));
hyp = packupScene(res.BESTHYPS(1), CAN_POOL{res.BESTROOM(1)});
%%
load(fullfile(config.bufname, config.localResultFile));
for i = 1:6
    for j = 1:3
        score(j) = ITERHYPS{j,i}.Scrs(1);
    end
    [~,bestCase] = max(score);
    AllScene(i) = packupScene(ITERHYPS{bestCase,i}.Hyps(1), CAN_POOL{ITERHYPS{bestCase,i}.roomid});
end

rot_file_name = [config.bufname config.IMGLIST(aid).name config.rotateImageSuffix];
rotImg = imread(rot_file_name);

figure;
for i = 1:6
    subplot(2,3,i);
    imshow(visualObject2D(rotImg, AllScene{i}));
end