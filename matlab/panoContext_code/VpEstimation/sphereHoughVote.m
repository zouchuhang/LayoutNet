function [ refiXYZ, lastStepCost, lastStepAngle ] = sphereHoughVote( segNormal, segLength, segScores, binRadius, orthTolerance, candiSet, force_unempty )
%SPHEREHOUGHVOTE Summary of this function goes here
%   Detailed explanation goes here
if ~exist('force_unempty','var')
    force_unempty = true;
end

%% initial guess
% segNormal = arcList(:,1:3);
% segLength = sqrt( sum((arcList(:,4:6)-arcList(:,7:9)).^2, 2));
% segScores = arcList(:,end);
numLinesg = size(segNormal,1);

% [voteBinPoints tri] = icosahedron2sphere(level);
voteBinPoints = candiSet;
voteBinPoints(voteBinPoints(:,3)<0,:) = []; 
reversValid = segNormal(:,3)<0;
segNormal(reversValid,:) = -segNormal(reversValid,:);

voteBinUV = xyz2uvN(voteBinPoints);
numVoteBin = length(voteBinPoints);
voteBinValues = zeros(numVoteBin,1);
for i = 1:numLinesg
    tempNorm = segNormal(i,:);
    tempDots = dot(voteBinPoints, repmat(tempNorm, [numVoteBin 1]), 2);
    
%     tempAngs = acos(abs(tempDots));
%     voteBinValues = voteBinValues + normpdf(tempAngs, 0, 0.5*binRadius*pi/180)*segScores(i)*segLength(i);
%     voteBinValues = voteBinValues + max(0, (2*binRadius*pi/180-tempAngs)./(2*binRadius*pi/180))*segScores(i)*segLength(i);
    
    
    valid = abs(tempDots)<cos((90-binRadius)*pi/180);

    voteBinValues(valid) = voteBinValues(valid) + segScores(i)*segLength(i);

end

checkIDs1 = find(voteBinUV(:,2)>pi/3);
voteMax = 0;
checkID1Max = 0;
checkID2Max = 0;
checkID3Max = 0;

for j = 1:length(checkIDs1)
%     fprintf('%d/%d\n', j, length(checkIDs1));
    checkID1 = checkIDs1(j); 
    vote1 = voteBinValues(checkID1);
    if voteBinValues(checkID1)==0 && force_unempty
        continue;
    end
    checkNormal = voteBinPoints(checkID1,:);
    dotProduct = dot(voteBinPoints, repmat(checkNormal, [size(voteBinPoints,1) 1]), 2);
    checkIDs2 = find(abs(dotProduct)<cos((90-orthTolerance)*pi/180));
    for i = 1:length(checkIDs2)
        checkID2 = checkIDs2(i);
        if voteBinValues(checkID2)==0 && force_unempty
            continue;
        end
        vote2 = vote1 + voteBinValues(checkID2);
        cpv = cross(voteBinPoints(checkID1,:), voteBinPoints(checkID2,:));
        cpn = sqrt(sum(cpv.^2));
        dotProduct = dot(voteBinPoints, repmat(cpv, [size(voteBinPoints,1) 1]), 2)./cpn;
        checkIDs3 = find(abs(dotProduct)>cos(orthTolerance*pi/180));    
        for k = 1:length(checkIDs3)
            checkID3 = checkIDs3(k); 
            if voteBinValues(checkID3)==0 && force_unempty
                continue;
            end
            vote3 = vote2 + voteBinValues(checkID3);
            if vote3>voteMax
%                 fprintf('%f\n', vote3);
                lastStepCost = vote3-voteMax;
                if voteMax ~= 0
                    tmp = dot(voteBinPoints([checkID1Max checkID2Max checkID3Max],:), ...
                              voteBinPoints([checkID1 checkID2 checkID3],:), 2);
                    lastStepAngle = acos(tmp);
                else
                    lastStepAngle = [0 0 0];
                end
                                     
                checkID1Max = checkID1;
                checkID2Max = checkID2;
                checkID3Max = checkID3;               
                                           
                voteMax = vote3;
            end
%             voteBins(checkID1, checkID2, checkID3) = true;
        end
    end
end

if checkID1Max==0
    fprintf('Warning: No orthogonal voting exist!!!\n');
    refiXYZ = [];
    lastStepCost = 0;
    lastStepAngle = 0;
    return;
end
initXYZ = voteBinPoints([checkID1Max checkID2Max checkID3Max],:);

%% refine
% binRadius = binRadius/2;

refiXYZ = zeros(3,3);
dotprod = dot(segNormal, repmat(initXYZ(1,:), [size(segNormal,1) 1]), 2);
valid = abs(dotprod)<cos((90-binRadius)*pi/180);
validNm = segNormal(valid,:);
validWt = segLength(valid).*segScores(valid);
validWt = validWt./max(validWt);
[~,refiNM] = curveFitting(validNm, validWt);
refiXYZ(1,:) = refiNM;

dotprod = dot(segNormal, repmat(initXYZ(2,:), [size(segNormal,1) 1]), 2);
valid = abs(dotprod)<cos((90-binRadius)*pi/180);
validNm = segNormal(valid,:);
validWt = segLength(valid).*segScores(valid);
validWt = validWt./max(validWt);
validNm(end+1,:) = refiXYZ(1,:);
validWt(end+1,:) = sum(validWt)*0.1;
[~,refiNM] = curveFitting(validNm, validWt);
refiXYZ(2,:) = refiNM;

refiNM = cross(refiXYZ(1,:), refiXYZ(2,:), 2);
refiXYZ(3,:) = refiNM./norm(refiNM);



% [~,refiNM] = curveFitting(validNm, validWt);
% refiXYZ(i,:) = refiNM;
% 
% 
% 
% for i = 1:3
%     dotprod = dot(segNormal, repmat(initXYZ(i,:), [size(segNormal,1) 1]), 2);
%     valid = abs(dotprod)<cos((90-binRadius)*pi/180);
%     validNm = segNormal(valid,:);
%     validWt = segLength(valid).*segScores(valid);
%     [~,refiNM] = curveFitting(validNm, validWt);
%     refiXYZ(i,:) = refiNM;
% end
%% output 
% % [voteBinPoints tri] = icosahedron2sphere(level);
% OBJ.vertices = voteBinPoints;
% % OBJ.vertices_normal = zeros(size(voteBinPoints));
% % OBJ.vertices_normal(:,1) = 1;
% OBJ.vertices_normal = voteBinPoints;
% uv = xyz2uvN(voteBinPoints);
% OBJ.vertices_texture = [(uv(:,1)+pi)/2/pi (uv(:,2)+pi/2)/pi];
% 
% OBJ.objects.type = 'f';
% OBJ.objects.data.vertices = tri;
% OBJ.objects.data.texture = tri;
% OBJ.objects.data.normal = tri;
% 
% % check boundary
% newVTSID = size(OBJ.vertices_texture,1);
% newFace = 0;
% newAddVT = zeros(1000,2);
% for i = 1:size(OBJ.objects.data.vertices,1)
%     texture = OBJ.objects.data.texture(i,:);
%     vt = OBJ.vertices_texture(texture,:);
%     v = OBJ.vertices(texture,:);
%     if (std(vt(:,1)))<0.3
%         continue;
%     end
%     
%     newFace = newFace + 1;
%     
%     modify = (vt(1,1)-vt(:,1))>0.5;
%     vt(modify,1) = vt(modify,1)+1;
%     modify = (vt(1,1)-vt(:,1))<-0.5;
%     vt(modify,1) = vt(modify,1)-1;
%     
%     newAddVT((newFace-1)*3+1:newFace*3,:) = vt;
%     OBJ.objects.data.texture(i,:) = [(newFace-1)*3+1:newFace*3] + newVTSID;
%     
%     if newFace>300
%         fprintf('Warning: pre-assign more memory!/n');
%     end
% end
% OBJ.vertices_texture = [OBJ.vertices_texture;newAddVT];
% 
% material(1).type='newmtl';
% material(1).data='Textured';
% material(2).type='Ka';
% material(2).data=[1.0 1.0 1.0];
% material(3).type='Kd';
% material(3).data=[1.0 1.0 1.0];
% material(4).type='Ks';
% material(4).data=[1.0 1.0 1.0];
% material(5).type='illum';
% material(5).data=2;
% material(6).type='map_Kd';
% material(6).data='pano_hotel_2.jpg';
% OBJ.material = material;
% 
% write_wobj(OBJ, 'exciting_notext.obj');

% uv = xyz2uvN(voteBinPoints);
% textureMap = zeros(512, 1024);
% x = min(round((uv(:,1)+pi)/2/pi*1024+1), 1024);
% y = 512 - min(round((uv(:,2)+pi/2)/pi*512+1), 512) + 1;
% value = voteBinValues./max(voteBinValues);
% value = value.^3;
% % [gridX, gridY] = meshgrid(x,y);
% for i = 1:length(x)
%     sx = max(1, x(i)-3);
%     ex = min(1024, x(i)+3);
%     sy = max(1, y(i)-3);
%     ey = min(512, y(i)+3);
%     textureMap(sy:ey, sx:ex) = value(i);
% end
% 
% % ind = sub2ind([512 1024], y, x);
% % textureMap(ind) = voteBinValues;
% imwrite(textureMap, 'exciting.png');

end

