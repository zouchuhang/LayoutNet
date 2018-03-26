function [ highPixelBox, highFeatureBox, feature] = rectDetectionFlex( img, model, config)
%EVALUATERECTDETECTOR Summary of this function goes here
%   Detailed explanation goes here
pattern = config.pattern;
patternSz = config.patternSz;
nmsThresh = config.nmsThresh;
repThresh = config.repThresh;
max_scale = config.max_scale;
min_scale = config.min_scale;


partTemplate = reshape(model.belta, 5, 5, 31, 9);
patternNum = size(pattern,1);

% for i = 1:patternNum
%     [template{i} tMask{i}] = feature2pattern( model.belta, pattern(i,:), patternSz(i,:), 5, 9);
% end

img = double(img);
% compute response
feature = feature_pyramid( img, max_scale, min_scale );
for pid = 1:9
    for sid = 1:length(feature.feat)
        match{pid, sid} = fconvblas(feature.feat{sid}, {partTemplate(:,:,:,pid)}, 1, 1);
        match{pid, sid} = match{pid, sid}{1};
    end
end
% combine score
for pid = 1:patternNum
    for sid = 1:length(feature.feat)
        [feaH, feaW, ~] = size(feature.feat{sid});
        resH = feaH-patternSz(pid,2)+1;
        resW = feaW-patternSz(pid,1)+1;
        response = zeros( resH, resW);
        for i = 1:9
            px = pattern(pid,2*i-1);
            py = pattern(pid,2*i);
            response = response + match{i,sid}(py:py+resH-1, px:px+resW-1);
        end
        report{pid,sid} = response;
    end
end
% get high report
for pid = 1:patternNum
    for sid = 1:length(feature.feat)
        bbs{pid,sid} = zeros(0, 5);
        if isempty(report{pid, sid})
            continue;
        end
        [scores,indexes] = sort(report{pid,sid}(:),'descend');
        if numel(scores)>5000
            choosen = 1:5000;
            indexes = indexes(choosen);
            scores  = scores(choosen);
        end
        [indexes1 indexes2] = ind2sub(size(report{pid,sid}),indexes);
        bbs{pid,sid} = ([indexes2-1+2.5 indexes1-1+2.5 indexes2+patternSz(pid,1)-1-2.5 indexes1+patternSz(pid,2)-1-2.5]) * 8 / feature.scale(sid) + repmat([1 1 0 0],numel(indexes),1);
        bbs{pid,sid} = [bbs{pid,sid} scores];
        fbs{pid,sid} = [indexes2 indexes1 indexes2+patternSz(pid,1)-1 indexes1+patternSz(pid,2)-1 sid*ones(length(indexes1),1) pid*ones(length(indexes1),1)];            
        
%         ratio = rectIntersection([344 836 774 1196], bbs{pid,sid});
%         [S, I] = max(ratio);
%         if S>0.8
%             fprintf('%d: %f %f %f %f %f\n', pid, bbs{pid,sid}(I,1),bbs{pid,sid}(I,2),bbs{pid,sid}(I,3),bbs{pid,sid}(I,4),bbs{pid,sid}(I,5));
%         end
    end
end
% get nms response
tic;
nmsBBS = zeros(0,5);
nmsFBS = zeros(0,6);
for pid = 1:patternNum
    scaleBBS = zeros(0,5);
    scaleFBS = zeros(0,6);
    for sid = 1:length(feature.feat)
        if isempty(bbs{pid, sid})
            continue;
        end
        valid = bbs{pid,sid}(:,5)>repThresh + model.b;
        tempBBS = bbs{pid,sid}(valid,:);
        tempFBS = fbs{pid,sid}(valid,:);
        indexes = nmsMe(tempBBS, nmsThresh);
        scaleBBS = [scaleBBS; tempBBS(indexes,:)];
        scaleFBS = [scaleFBS; tempFBS(indexes,:)];
    end
    indexes = nmsMe(scaleBBS, nmsThresh);
    nmsBBS = [nmsBBS; scaleBBS(indexes,:)];
    nmsFBS = [nmsFBS; scaleFBS(indexes,:)];
end
toc;

tic;
indexes = nmsMe(nmsBBS, nmsThresh);
toc;

highPixelBox = nmsBBS(indexes,:);
highFeatureBox = nmsFBS(indexes,:);

highPixelBox(:,5) = highPixelBox(:,5) - model.b;
valid = highPixelBox(:,5) >= repThresh;
highPixelBox = highPixelBox(valid, :);
highFeatureBox = highFeatureBox(valid, :);
% figure; imshow(img./255); hold on
% color = 'gcybrrrr';
% for i = 1:min(8, size(highPixelBox,1))
%     b = highPixelBox(i,:);
%     line(b([1 3 3 1 1]), b([2 2 4 4 2]), 'color', color(i), 'LineWidth', 2);
%     text(b(1),b(2)-15,sprintf('%f',b(5)),'color',color(i));
% end
% % visualize
% figure;
% for i = 1:size(highPixelBox,1)
%     f = highFeatureBox(i,:);
%     subplot(2,2,1); showHOG(template{f(6)});
% %     subplot(2,2,2); showHOG(feature.feat{5});
%     subplot(2,2,3); showHOG(feature.feat{f(5)}(f(2):f(4), f(1):f(3), :));
%     subplot(2,2,4); imshow(img./255);
%     b = highPixelBox(i,:);
%     line(b([1 3 3 1 1]), b([2 2 4 4 2]), 'color', 'r', 'LineWidth', 2);
%     title(sprintf('Score: %f', b(5)));
%     pause;
% end

end

