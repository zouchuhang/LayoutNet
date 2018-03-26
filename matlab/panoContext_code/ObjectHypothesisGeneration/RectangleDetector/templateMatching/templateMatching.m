function [ highPixelBox, highFeatureBox ] = templateMatching( feature, template )
%TEMPLATEMATCHING Summary of this function goes here
%   Detailed explanation goes here
highPixelBox = [];
highFeatureBox = [];

for level=1:length(feature.feat)
    % convolution
    match{level} = fconvblas(feature.feat{level}, {template}, 1, 1);
    match{level} = match{level}{1};
    
    % pick the top 5000 results
    [scores,indexes] = sort(match{level}(:),'descend');
    if numel(scores)>5000
        choosen = 1:5000;
        indexes = indexes(choosen);
        scores  = scores(choosen);
    end
    
    % obtain window pixel locations
    [indexes1 indexes2] = ind2sub(size(match{level}),indexes);
    bbs{level} = ([indexes2 indexes1 indexes2+size(template,2) indexes1+size(template,1)] - 1) * 8 / feature.scale(level) + repmat([1 1 0 0],numel(indexes),1);
    bbs{level} = [bbs{level} scores];
    fbs{level} = [indexes2 indexes1 indexes2+size(template,2)-1 indexes1+size(template,1)-1 level*ones(length(indexes1),1)];
        
    % non maximal suppression to remove windows too close
    indexes = nmsMe(bbs{level}, 0.5); % pascal voc 0.5 criteria
    bbs{level} = bbs{level}(indexes,:);    
    fbs{level} = fbs{level}(indexes,:);
    highPixelBox = [highPixelBox; bbs{level}];
    highFeatureBox = [highFeatureBox; fbs{level}];
end

% overall best match
indexes = nmsMe(highPixelBox, 0.5); % pascal voc 0.5 criteria
highPixelBox = highPixelBox(indexes,:); 
highFeatureBox = highFeatureBox(indexes, :);

end

