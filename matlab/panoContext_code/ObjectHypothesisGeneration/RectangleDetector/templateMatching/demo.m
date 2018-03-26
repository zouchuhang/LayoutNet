% a simple demo of HOG template matching : Jianxiong 

% you may need to compile some mex files in matlab: see features_compile.m

clear
close all

image = imread('demo_image.jpg');

feature = feature_pyramid(image);

% pick a template
template_level = 1;
template_width = 8;
template_height= 6;
template_left  = 41;
tempalte_top   = 31;

% extract a tempalte
template = feature.feat{template_level}(tempalte_top:(tempalte_top+template_height-1),template_left:(template_left+template_width-1),:);

% normalize the tempalte to zero mean
% empirically better for some task, you can also try to disable this
template = template - mean(template(:));

% corresponding image
template_box = ([template_left tempalte_top template_left+size(template,2) tempalte_top+size(template,1)] - 1) * 8 / feature.scale(template_level) + [1 1 0 0];
template_box = round(template_box);
tempalte_image = image(template_box(2):template_box(4),template_box(1):template_box(3),:);

matchMatrix = [];

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
    
    % non maximal suppression to remove windows too close
    indexes = nmsMe(bbs{level}, 0.5); % pascal voc 0.5 criteria
    bbs{level} = bbs{level}(indexes,:);    
    
    matchMatrix = [matchMatrix; bbs{level}];
end

% overall best match
indexes = nmsMe(matchMatrix, 0.5); % pascal voc 0.5 criteria
matchMatrix = matchMatrix(indexes,:);  

%% visualization

visLevel = [];
for level = 1:3:length(feature.feat)
    if ~isempty(bbs{level})
        visLevel = [visLevel, level];
    end
end

subplot(length(visLevel)+1,3, 1);
showHOG(template);
title('the HOG tempalte');
subplot(length(visLevel)+1,3, 2);
imshow(tempalte_image);
title('image of the template');
subplot(length(visLevel)+1,3, 3);
imshow(image);
hold on
for b=1:min(5,size(matchMatrix))
    plot(matchMatrix(b,[1 3 3 1 1]),matchMatrix(b,[2 2 4 4 2]),'-y');
end
title('top 5 detection over all scale levels');

for l=1:length(visLevel)
    level = visLevel(l);
    subplot(length(visLevel)+1,3, 3*l+1);
    showHOG(feature.feat{level});
    title(sprintf('Level %d HOG ',level));
    
    subplot(length(visLevel)+1,3, 3*l+2);
    imagesc(match{level});
    axis image
    axis off
    title(sprintf('Level %d response ',level));
    
    subplot(length(visLevel)+1,3, 3*l+3);
    imshow(image);
    hold on
    for b=1:min(5,size(bbs{level},1)) % only visualization the top five results
        plot(bbs{level}(b,[1 3 3 1 1]),bbs{level}(b,[2 2 4 4 2]),'-y');
    end
    title(sprintf('Level %d top 5',level));
end
 

