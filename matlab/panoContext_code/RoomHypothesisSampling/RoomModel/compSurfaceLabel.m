function [ sepScene, wallPanoNormal] = compSurfaceLabel( rotImg )%compSurfaceLabel( data, folderName )
%COMPSURFACELABEL Compute geometric context surface label for panorama
%   sepScene: surface label for each seprate view, 6 channels
%   wallPanoNormal: surface label for panorama, 7 channels

% [rotImg, data, R] = loadBuffer( data, folderName ); 
%% emprical good separate of panorama
cutSize = 640;
% fov = pi*2/3;
fov = [pi/2*ones(1,8) pi/2*ones(1,8)];
xh = [-pi:(pi/4):(3/4*pi) -pi:(pi/4):(3/4*pi)];
yh = [zeros(1, 8) -1/6*pi*ones(1,8)];

% fov = [pi/2*ones(1,36) pi/2*ones(1,36)];
% xh = [-pi:(pi/18):(17/18*pi) -pi:(pi/18):(17/18*pi)];
% yh = [zeros(1,36) -1/6*pi*ones(1,36)];

% fov = [pi/2*ones(1,8)];
% xh = [-pi:(pi/4):(3/4*pi)];
% yh = [zeros(1,8)];

% fov = [pi/2*ones(1,36)];
% xh = [-pi:(pi/18):(17/18*pi)];
% yh = [zeros(1,36)];
[sepScene] = separatePano( rotImg, fov, xh, yh, cutSize);

%% 
% imdir = './roomModel/SpatialLayout_shrink/Images_resized/';
% workspcdir = './roomModel/SpatialLayout_shrink/tempworkspace/';
% imsegdir = './roomModel/SpatialLayout_shrink/Imsegs/';
% moveFolder = './roomModel/SpatialLayout_shrink/';
sepScene(1).cimages = [];
sepScene(1).maxCimage = [];
numScene = length(sepScene);
for sid = 1:numScene
    img = sepScene(sid).img;
    img = imresize(img, [640 640]);
    try
        [ cimages ] = getSurfaceLabelInitial( img );
    catch
        cimages{1} = 1/7*ones(cutSize,cutSize,7);
    end
%     [cimages, Polyg] = getSurfaceLabelLayout( sepScene(sid), data  )
    sepScene(sid).cimages = cimages;  
    [~, indd] = max(cimages{1}, [], 3);
    sepScene(sid).maxCimage = indd;
end

%% combine
% 1 middle wall=yellow; 2 left wall=red; 3 right wall=cyan; 4 floor=green; 5 ceiling=blue; 6 objects=pink
width = 2048;   height = 1024;
wallPano = zeros(height, width, 7);
wallPanoWeight = zeros(height, width);
for vid = 1:length(sepScene)
    vx = sepScene(vid).vx;
    
    [curPano, curPanoValid] = ...
    im2Sphere( sepScene(vid).cimages{1}, sepScene(vid).fov, width, height, ...
                sepScene(vid).vx, sepScene(vid).vy);
    curPano(~curPanoValid) = 0;
    wallPanoWeight(curPanoValid(:,:,1)) = wallPanoWeight(curPanoValid(:,:,1)) + 1;
    
    normMiddle = vx+pi;
    normMLR = normValue([normMiddle normMiddle-pi/2 normMiddle+pi/2], -pi, pi);
    normGlobal = [-pi -pi/2 0 pi/2];
    for i = 1:3
        angle = abs(normValue(normGlobal - normMLR(i), -pi, pi));
        angle(angle>pi/2) = pi/2; % cos^2 weight
%         [~,I] = sort(angle, 'ascend'); % nearest neighbor
%         angle(I(1)) = 0;
%         angle(I(2:4)) = pi/2;
        for j = find(angle<pi/2)
            wallPano(:,:,j) = wallPano(:,:,j) + curPano(:,:,i)*(cos(angle(j)).^2);
        end
    end
    
    wallPano(:,:,5) = wallPano(:,:,5) + curPano(:,:,4); % floor
    wallPano(:,:,6) = wallPano(:,:,6) + curPano(:,:,5); % ceiling
    wallPano(:,:,7) = wallPano(:,:,7) + curPano(:,:,6); % object
end
wallPanoNormal = wallPano./repmat(wallPanoWeight+0.000001,[1 1 7]);



%% more informative combine
% % Left wall: view+(90-FOV/2) ~ view+(90-FOV/2)+180
% % Right wall: view-(90-FOV/2)-180 ~ view-(90-FOV/2)
% % Middle wall: Left wall + Right wall
% width = 1024;   height = 512;
% wallPano = zeros(512, 1024, 4);
% for vid = 1:length(sepScene)
%     vx = sepScene(vid).vx;
%     fov = sepScene(vid).fov;
%     [curPano, curPanoValid] = ...
%         im2Sphere( sepScene(vid).cimages{1}, sepScene(vid).fov, width, height, ...
%                     sepScene(vid).vx, sepScene(vid).vy);
%     curPano(~curPanoValid) = 0;
%     %Left
%     lowBound = vx+(pi/2-fov/2);
%     uppBound = vx+(pi/2-fov/2)+pi;
%     range = normRange(lowBound, uppBound, -pi, pi);
%     Linside = insideRange([-pi -pi/2 0 pi/2], range);
% %     Lweight = [1 1 1 1];
%     %Right
%     lowBound = vx-(pi/2-fov/2)-pi;
%     uppBound = vx-(pi/2-fov/2);
%     range = normRange(lowBound, uppBound, -pi, pi);
%     Rinside = insideRange([-pi -pi/2 0 pi/2], range);
% %     Rweight = [1 1 1 1];
%     %Middle
%     lowBound = vx+(pi/2-fov/2);
%     uppBound = vx-(pi/2-fov/2);
%     range = normRange(lowBound, uppBound, -pi, pi);
%     Minside = insideRange([-pi -pi/2 0 pi/2], range);
% %     Mweight = [1 1 1 1];
%     %Vote Left
%     score = curPano(:,:,2);
%     voteMatrix = repmat(reshape(Linside, [1 1 4]), [512 1024 1]);
%     voteMatrix = voteMatrix .* repmat(score, [1 1 4]);
%     wallPano = wallPano + voteMatrix;
%     
%     %Vote Right
%     score = curPano(:,:,3);
%     voteMatrix = repmat(reshape(Linside, [1 1 4]), [512 1024 1]);
%     voteMatrix = voteMatrix .* repmat(score, [1 1 4]);
%     wallPano = wallPano + voteMatrix;
%     
%     %Vote Middle
%     score = curPano(:,:,1);
%     voteMatrix = repmat(reshape(Linside, [1 1 4]), [512 1024 1]);
%     voteMatrix = voteMatrix .* repmat(score, [1 1 4]);
%     wallPano = wallPano + voteMatrix;
%     
% end

%% heuristic combine
% MASK = zeros(640,640);
% MASK(213:427,1:320) = 1;
% MASK(213:427,321:640) = 2;
% All_MLR2ABCD = zeros(8,3);
% for vid = [2 4 6 8]
%     histLabelL = hist(sepScene(vid).maxCimage(MASK<1.5 & MASK>0.5), 1:3);
%     histLabelR = hist(sepScene(vid).maxCimage(MASK>1.5 & MASK<2.5), 1:3);
%     histLabelL = histLabelL./sum(histLabelL) + 0.0001;
%     histLabelR = histLabelR./sum(histLabelR) + 0.0001;
%     entropyL = -sum(histLabelL.*log2(histLabelL));
%     entropyR = -sum(histLabelR.*log2(histLabelR));
%     MLR2ABCD = [0 0 0]; %[M L R], 1234 for ABCD, 0 for Null
%     histLabel = hist(sepScene(vid).maxCimage(MASK<2.5 & MASK>0.5), 1:3);
%     [~, IW] = sort(histLabel, 'descend');
%     if entropyL<entropyR
%         [~, I] = max(histLabelL);
%         MLR2ABCD(I) = vid/2;
%         if ~ismember(I, IW(1:2))
%             fprintf('Probably some problems\n');
%         end
%         if I==1
%             MLR2ABCD(3) = rem(vid/2,4)+1;
%         else
%             MLR2ABCD(setdiff(IW(1:2),I)) = rem(vid/2,4)+1;
%         end
%     else
%         [~, I] = max(histLabelR);
%         MLR2ABCD(I) = rem(vid/2,4)+1;
%         if ~ismember(I, IW(1:2))
%             fprintf('Probably some problems\n');
%         end
%         if I==1
%             MLR2ABCD(2) = vid/2;
%         else
%             MLR2ABCD(setdiff(IW(1:2),I)) = vid/2;
%         end
%     end
%     All_MLR2ABCD(vid,:) = MLR2ABCD;
% end
% width = 1024;   height = 512;
% for vid = [1 3 5 7]
%     preID = vid - 1; preID(preID==0) = 8;
%     latID = vid + 1;
%     
%     curMaxCimage = sepScene(vid).maxCimage;
%     curMaxCimage(MASK<0.5) = 0;
%     [curPano, curPanoValid] = ...
%         im2Sphere( curMaxCimage, sepScene(vid).fov, width, height, ...
%                     sepScene(vid).vx, sepScene(vid).vy);
%     curPano(~curPanoValid) = 0;
%     
%     preMaxCimage = sepScene(preID).maxCimage;
%     preMaxCimage(MASK<0.5) = 0;
%     [prePano, prePanoValid] = ...
%         im2Sphere( preMaxCimage, sepScene(preID).fov, width, height, ...
%                     sepScene(preID).vx, sepScene(preID).vy);
%     prePano(~prePanoValid) = 0;    
%     
%     latMaxCimage = sepScene(latID).maxCimage;
%     latMaxCimage(MASK<0.5) = 0;
%     [latPano, latPanoValid] = ...
%         im2Sphere( latMaxCimage, sepScene(latID).fov, width, height, ...
%                     sepScene(latID).vx, sepScene(latID).vy);
%     latPano(~latPanoValid) = 0;    
%     
%     preConsist = [0 0 0];
%     preNumber = [0 0 0];
%     for i = 1:3
%         mask = curPano>i-0.5 & curPano<i+0.5 & prePano>0.5 & prePano<3.5;
%         h = hist(prePano(mask),1:3);
%         [B,I] = max(h);
%         J = All_MLR2ABCD(preID,I);
%         preConsist(i) = J;
%         preNumber(i) = B;
%         
%         mask = curPano>i-0.5 & curPano<i+0.5;
%         if sum(mask(:))<250 || B/sum(mask(:))<0.05
%             preConsist(i) = 0;
%         end
%     end
%     latConsist = [0 0 0];
%     latNumber = [0 0 0];
%     for i = 1:3
%         mask = curPano>i-0.5 & curPano<i+0.5 & latPano>0.5 & latPano<3.5;
%         h = hist(latPano(mask),1:3);
%         [B,I] = max(h);
%         J = All_MLR2ABCD(latID,I);
%         latConsist(i) = J;
%         latNumber(i) = B;
%         
%         mask = curPano>i-0.5 & curPano<i+0.5;
%         if sum(mask(:))<250 || B/sum(mask(:))<0.05
%             latConsist(i) = 0;
%         end
%     end
%     for i = 1:3
%         if preNumber(i)>latNumber(i)
%             All_MLR2ABCD(vid, i) = preConsist(i);
%         else
%             All_MLR2ABCD(vid, i) = latConsist(i);
%         end
%     end
% end
% wallPano = zeros(512, 1024, 4);
% for vid = [2 4 6 8]
%     [curPano, curPanoValid] = ...
%     im2Sphere( sepScene(vid).cimages{1}, sepScene(vid).fov, width, height, ...
%                 sepScene(vid).vx, sepScene(vid).vy);
%     curPano(~curPanoValid) = 0;
%     voteMat = zeros(512, 1024, 4);
%     for i = find(All_MLR2ABCD(vid,:)~=0)
%         voteMat(:,:,All_MLR2ABCD(vid,i)) = curPano(:,:,i);
%     end
%     wallPano = wallPano + voteMat;
% end


end

function [nv] = normValue( v, minBound, maxBound)
nv = rem((v-minBound)+2*(maxBound-minBound), maxBound-minBound) + minBound;
end

function [range] = normRange( low, upp, minBound, maxBound)
normLow = rem((low-minBound)+2*(maxBound-minBound), maxBound-minBound) + minBound;
normUpp = rem((upp-minBound)+2*(maxBound-minBound), maxBound-minBound) + minBound;
if normLow<=normUpp
    range = [normLow normUpp];
else
    range = [normLow maxBound; minBound normUpp];
end
end

function [inside] = insideRange(values, range)
inside = false(size(values));
for i = 1:size(range,1)
    inside = inside | (values>=range(i,1) & values<=range(i,2));
end
end
