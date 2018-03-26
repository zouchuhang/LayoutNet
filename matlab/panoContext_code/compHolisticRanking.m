function result = compHolisticRanking( aid )
%COMPHOLISTICRANKING Summary of this function goes here
%   Detailed explanation goes here
global config
bufname = config.bufname;

did = config.valid_data_id(config.valid_anno_id==aid);

load(config.globalSceneSVMFile);
scenesvm = MODEL{config.globalSceneSVMiterID}.initModel;
clear MODEL

load(sprintf([bufname config.sceneFeatFile], 1));
SEED_HYPS = ALL_HYPS(1);
names = fieldnames(SEED_HYPS);
for i = 1:length(names)
    SEED_HYPS.(names{i}) = [];
end

HYPS = repmat(SEED_HYPS, 5000, 1);
ALL_COST = -10*ones(5000,1);
for rid = 1:config.sampleroomnum
    fprintf('>>%d:%d ', aid, rid);
    try
        load(sprintf([bufname '/' config.sceneFeatFile], rid));
    catch
        fprintf('no data\n');
        continue;
    end
    if isempty(ALL_HYPS)
        continue;
    end
    HYPS((rid-1)*100+1:rid*100) = ALL_HYPS;
    COST((rid-1)*100+1:rid*100) = vertcat(ALL_HYPS.COST);


    [ ~, hypScore] = sceneHypsEvaluation( vertcat(ALL_HYPS.sceneImgFea), ...
            vertcat(ALL_HYPS.MINSCORE), scenesvm, did );
    ALL_COST((rid-1)*100+1:rid*100) = hypScore;  
    ROOMID((rid-1)*100+1:rid*100) = rid;
end
fprintf('\n');

load([bufname config.objHypoFile3D]);

%% get our detection performance
[~,I] = sort(ALL_COST,'descend');
%keyboard
scene = packupScene(HYPS(I(1)), CAN_POOL{ROOMID(I(1))});
DET = scene{1};

result.Hyps = HYPS(I(1));
result.GLOBAL = ALL_COST(I(1));
result.ROOMID = ROOMID(I(1));
result.DET = DET;

rot_file_name = [bufname config.IMGLIST(aid).name config.rotateImageSuffix];
rotImg = imread(rot_file_name);
figure; imshow(visualObject2D( rotImg, DET, 512, 1024 ));
figure; visualObject3D( DET, 4, 1 );

end

