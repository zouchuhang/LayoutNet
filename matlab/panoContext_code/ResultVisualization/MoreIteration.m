BEDROOMDATA;

ResultDir = '/n/fs/sun360/panoParser/SceneParsing/topdowninference/scenefeature/localsampleresult/';
% load(fullfile(config.folderName, config.objHypoFile3D));

OutputDir = '../MoreIteration/';
if ~exist(OutputDir, 'dir')
    mkdir(OutputDir);
end

% write HTML
% HTMLDir = fullfile(OutputDir, 'HTML');
HTMLDir = OutputDir;
if ~exist(HTMLDir, 'dir')
    mkdir(HTMLDir);
end
webpagename = fullfile(HTMLDir,'longiter.html');
copyfile('/n/fs/sun360/web2014/LocalSampling/tabletemplate.html', webpagename);
FID = fopen(webpagename,'a');
fprintf(FID, '<tr>\n');

for did = 1:length(config.valid_anno_id)
    aid = config.valid_anno_id(did);
    fprintf('did: %d, aid: %d\n', did, aid);

    load(fullfile(config.folderName, config.IMGLIST(aid).name, config.objHypoFile3D));
    
    rotImg = imread(fullfile(config.folderName, config.IMGLIST(aid).name, ...
        [config.IMGLIST(aid).name config.rotateImageSuffix]));
    
    if aid<=87
        load(sprintf([ResultDir '%03d_dense_new_v2_global_iter6.mat'], aid));
    else
        load(sprintf([ResultDir '%03d_dense_new_v2_global_iter3.mat'], aid));
    end
    [hypNum, iterNum] = size(ITERHYPS);
    
%     for iter = 1:iterNum
%         for hid = 1:hypNum
%             if isempty(ITERHYPS{hid,iter}) || isempty(ITERHYPS{hid,iter}.Hyps)
%                 continue;
%             end
%             ITERSCENE = packupScene(ITERHYPS{hid,iter}.Hyps(1), CAN_POOL{ITERHYPS{hid,iter}.roomid});
%             visual2D = visualObject2D(rotImg, ITERSCENE{1});
%             imwrite(visual2D, fullfile(OutputDir, sprintf('best_%d_%d_%d.jpg',aid,hid,iter)));
%         end
%     end
 
    imgpos = fullfile('./data', config.IMGLIST(aid).name, ...
        [config.IMGLIST(aid).name config.rotateImageSuffix]);

    fprintf(FID, '\t<td rowspan="%d">%03d</td>\n', hypNum, aid);
    fprintf(FID, '\t<td rowspan="%d"><img class="s1" src="%s" alt="" width="512" height="256">%s</td>\n', ...
        hypNum, imgpos, config.IMGLIST(aid).name);
    
    for hid = 1:hypNum
        for iter = 1:iterNum
            if isempty(ITERHYPS{hid,iter}) || isempty(ITERHYPS{hid,iter}.Hyps)
                continue;
            end
            fprintf(FID, '<td><img class="s1" src="./best_%d_%d_%d.jpg" alt="" width="512" height="256">%f</td>\n', ...
                aid, hid, iter, ITERHYPS{hid, iter}.Scrs(1));
        end
        fprintf(FID, '</tr>\n');
        fprintf(FID, '<tr>\n');
    end
    
end