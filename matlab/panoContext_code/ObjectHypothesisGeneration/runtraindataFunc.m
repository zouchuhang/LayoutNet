function runtraindataFunc( aid )
%RUNTRAINDATAFUNC Summary of this function goes here
%   Detailed explanation goes here
% matlabpool(2);
% matlabpool open local 4
% aid = aid+83;

% IDS = [7,10,22,28,31,33,35,36,37,42,71,72,75,77,81,89,91,92,101,112,114,122,124,132,138,156,157,158,160,162,168,174,175,187,193,201,204,207,223,226,230,237,250,252,255,259,270,271,274,276,289,299,308,317,341,342,347,350,355,362,363,372,381,382,383,389,391,392,401,423];
% IDS = [184];
% aid = IDS(aid);

[~,b] = system('hostname');
fprintf('Running on %s\n',b);
t = clock();
fprintf('Run at %d.%d.%d, %d:%d\n', t(1:5));

fprintf('aid: %d\n', aid);

addpath(genpath('/n/fs/sun360/panoParser/SceneParsing'));
cd '/n/fs/sun360/panoParser/SceneParsing';

folderName = './annotation/dataset/bedroom/';
load([folderName '/IMGLIST.mat']);
load([folderName '/ANNO_ALL.mat']);
load([folderName '/423GndTransform.mat']);

% testfile2 = [folderName IMGLIST(aid).name '/RECALL_morescore.mat'];
% if exist(testfile2, 'file')
% 	fprintf('Job: %d: %s\nPreviously done^_^\n', aid, testfile2);
% 	return;
% end

anno_valid = false( 1, length(ANNO_ALL));
for i = 1:length(anno_valid)
    anno_valid(i) = ~isempty(ANNO_ALL(i).ANNO3D) && ANNO_ALL(i).ANNO3D.b_singleroom;
end
if ~anno_valid(aid)
    return;
end

% valid_gnds = true(length(TRANS_GNDS),1);
% for i = 1:length(TRANS_GNDS)
%     if isempty(TRANS_GNDS(i).SAMPLE_COMP)
%         valid_gnds(i) = false;
%     end
% end
% VALID_TRANS_GNDS = TRANS_GNDS(valid_gnds);
% fprintf('Processing training data: %d/%d\n', aid, numData);
% if ~valid_gnds(aid)
%     return;
% end


% load([folderName '/ALL_GNDS_ROT.mat']);
% [ ALL_GNDS_ROT, ~, ~, ~] = gndRoomProcess( folderName );
% fprintf('Transform GND Rooms\n');

% numData = length(ANNO_ALL);


% % temporary job: update old hypotheis to new feature
% load([folderName ANNO_ALL(aid).name '/HYPFEATURE_NEW.mat']);
% FEATURES = compWholeHypFeatureA( aid, ALL_HYPS, [], folderName, ANNO_ALL );
% save([folderName ANNO_ALL(aid).name '/HYPFEATURE_NEW_UPDATE.mat'],'ALL_HYPS', 'FEATURES', 'COST');

% numhyp = 50;
% numsample = 200;
% [ ALL_HYPS, CAN_POOL, HYP_ROOM_ID ] = generateSeeds_P1( aid, numhyp, numsample, IMGLIST, folderName );
% [GOOD_HYPS, ~] = sampleFitGT( CAN_POOL, ALL_GNDS_ROT{aid,1} );
% ALL_HYPS = [ALL_HYPS;GOOD_HYPS];
% fprintf('Generated Room Hypothesis\n');
CAN_POOL = generateObjHyps( aid, 20000, IMGLIST, folderName  );
save([folderName ANNO_ALL(aid).name '/RECALL_final_P1_more.mat'],'CAN_POOL', '-v7.3');
fprintf('finish.\n');

% I = cellfun(@isempty,ALL_HYPS);
% ALL_HYPS = ALL_HYPS(~I);
% HYP_COMP = getHypInfo(ALL_HYPS);
% FEATURES = compWholeHypFeatureA( aid, ALL_HYPS, [], folderName, ANNO_ALL );
% fprintf('Computed feature\n');    
% %     load([folderName ANNO_ALL(aid).name '/HYPFEATURE_NEW.mat']);
% %     FEATURES = compWholeHypFeature( aid, ALL_HYPS, [], folderName, ANNO_ALL );
% 
% load([folderName ANNO_ALL(aid).name '/vanishing_point.mat']);
% hyp_vanishing = vp(3:-1:1,:);
% gnd_vanishing = ANNO_ALL(aid).vp(3:-1:1,:);
% gnd = ALL_GNDS_ROT{aid,1};
% COST = 50*ones(length(ALL_HYPS),1);
% I2H = diag([1 1 1])/hyp_vanishing';
% I2G = ANNO_ALL(aid).ANNO3D.R * ANNO_ALL(aid).ANNO3D.Rc;
% H2G_R = I2G/I2H;
% for i = 1:length(ALL_HYPS)
%     fprintf('HYPID: %d\n',i);
%     if isempty(ALL_HYPS{i})
%         continue;
%     end
%     COST(i) = roomLossFunction3D( gnd, ALL_HYPS{i}, H2G_R );
% end 
% 
% MINSCORE = zeros(length(ALL_HYPS), length(VALID_TRANS_GNDS));
% MINTRANS = zeros(length(ALL_HYPS), length(VALID_TRANS_GNDS));
% tic;
% for j = 1:length(VALID_TRANS_GNDS)
%     fprintf('GND: %d/%d\n', j, length(VALID_TRANS_GNDS));
%     tCOST = roomAlignmentMex(HYP_COMP, VALID_TRANS_GNDS(j).SAMPLE_COMP, 3, 0);
%     [B,I] = min( tCOST, [], 2);
%     MINSCORE(:,j) = B;
%     MINTRANS(:,j) = I;
% end
% toc;
% fprintf('Finish\n');
% % best ever
% 
%     
% save([folderName ANNO_ALL(aid).name '/HYPFEATURE_NEW_P2.mat'],'ALL_HYPS', 'FEATURES', 'COST', 'HYP_ROOM_ID', 'numhyp', 'numsample');
% save([folderName ANNO_ALL(aid).name '/ROOM_MATCHING_P2.mat'], 'MINSCORE', 'MINTRANS');
% save([folderName ANNO_ALL(aid).name '/CANDIDATEPOOL_P2.mat'], 'CAN_POOL');


end

