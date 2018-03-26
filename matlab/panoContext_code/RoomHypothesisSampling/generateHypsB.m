function [ hyps ] = generateHypsB( lines, vp, tol, omap, gc )
%GENERATEHYPS generate hypothesis, screen out inconsistent ones according
%to OM, GC, and merged map
%   lines: line segment
%   vp: vanishing point, in z,y,x order
%   tol: tolerance when assigning lines to vp
%   omap & gc: evidence
%   hyps: all hyps with >0.5 consistency with OM will be returned;

%% >> ASSIGN LINE SEGMENTS WITH TYPE

%% get line segments align with vanishing point
lines1 = getVanishingLine(lines, vp(1,:), tol);
lines2 = getVanishingLine(lines, vp(2,:), tol);
lines3 = getVanishingLine(lines, vp(3,:), tol);

%% extend line representation with line type
% line: [1x3 n] [planeID] [1x2 uv1] [u/d f/b r/l] (1/-1, 0 not sure, 2 not available)
lines1_attr = zeros(size(lines1,1),9);
lines2_attr = zeros(size(lines2,1),9);
lines3_attr = zeros(size(lines3,1),9);
lines1_attr(:,1:6) = lines1(:,1:6);
lines2_attr(:,1:6) = lines2(:,1:6);
lines3_attr(:,1:6) = lines3(:,1:6);

%% horizontal line, line connecting vp2 and vp3, u=1, d=-1
horizontal = cross(vp(2,:),vp(3,:));
horizontal = horizontal./norm(horizontal);
if horizontal(:,3)<0
    horizontal = -horizontal;
end
% line belong to vp2 and vp3 should have this attribute
% lines1 not available
lines1_attr(:,7) = 2;
% lines2
for i = 1:size(lines2_attr,1)
    l = lines2(i,:);
    u = l(5:6)*2*pi-pi;
    v = computeUVN( l(1:3), u, l(4));
    xyz = uv2xyzN([u' v'], l(4));
    b1 = dot(xyz(1,:), horizontal)>0;
    b2 = dot(xyz(2,:), horizontal)>0;
    if b1&&b2
        lines2_attr(i,7) = 1;
    elseif ~b1&&~b2
        lines2_attr(i,7) = -1;
    else
        lines2_attr(i,7) = 0;
    end   
end
% lines3
for i = 1:size(lines3_attr,1)
    l = lines3(i,:);
    u = l(5:6)*2*pi-pi;
    v = computeUVN( l(1:3), u, l(4));
    xyz = uv2xyzN([u' v'], l(4));
    b1 = dot(xyz(1,:), horizontal)>0;
    b2 = dot(xyz(2,:), horizontal)>0;
    if b1&&b2
        lines3_attr(i,7) = 1;
    elseif ~b1&&~b2
        lines3_attr(i,7) = -1;
    else
        lines3_attr(i,7) = 0;
    end   
end

%% front and back, line connecting vp1 and vp3, f=1, b=-1
frontNback = cross(vp(1,:),vp(3,:));
frontNback = frontNback./norm(frontNback);
if frontNback(:,2)<0
    frontNback = -frontNback;
end
% line belong to vp1 and vp3 should have this attribute
% lines2 not available
lines2_attr(:,8) = 2;
% lines1
for i = 1:size(lines1_attr,1)
    l = lines1(i,:);
    u = l(5:6)*2*pi-pi;
    v = computeUVN( l(1:3), u, l(4));
    xyz = uv2xyzN( [u' v'], l(4));
    b1 = dot(xyz(1,:), frontNback)>0;
    b2 = dot(xyz(2,:), frontNback)>0;
    if b1&&b2
        lines1_attr(i,8) = 1;
    elseif ~b1&&~b2
        lines1_attr(i,8) = -1;
    else
        lines1_attr(i,8) = 0;
    end   
end
% lines3
for i = 1:size(lines3_attr,1)
    l = lines3(i,:);
    u = l(5:6)*2*pi-pi;
    v = computeUVN( l(1:3), u, l(4));
    xyz = uv2xyzN( [u' v'], l(4));
    b1 = dot(xyz(1,:), frontNback)>0;
    b2 = dot(xyz(2,:), frontNback)>0;
    if b1&&b2
        lines3_attr(i,8) = 1;
    elseif ~b1&&~b2
        lines3_attr(i,8) = -1;
    else
        lines3_attr(i,8) = 0;
    end   
end

%% left and right, line connecting vp1 and vp2, r=1, l=-1
leftNright = cross(vp(1,:),vp(2,:));
leftNright = leftNright./norm(leftNright);
if leftNright(:,1)<0
    leftNright = -leftNright;
end
% line belong to vp1 and vp2 should have this attribute
% lines3 not available
lines3_attr(:,9) = 2;
% lines1
for i = 1:size(lines1_attr,1)
    l = lines1(i,:);
    u = l(5:6)*2*pi-pi;
    v = computeUVN( l(1:3), u, l(4));
    xyz = uv2xyzN( [u' v'], l(4));
    b1 = dot(xyz(1,:), leftNright)>0;
    b2 = dot(xyz(2,:), leftNright)>0;
    if b1&&b2
        lines1_attr(i,9) = 1;
    elseif ~b1&&~b2
        lines1_attr(i,9) = -1;
    else
        lines1_attr(i,9) = 0;
    end   
end
% lines2
for i = 1:size(lines2_attr,1)
    l = lines2(i,:);
    u = l(5:6)*2*pi-pi;
    v = computeUVN( l(1:3), u, l(4));
    xyz = uv2xyzN( [u' v'], l(4));
    b1 = dot(xyz(1,:), leftNright)>0;
    b2 = dot(xyz(2,:), leftNright)>0;
    if b1&&b2
        lines2_attr(i,9) = 1;
    elseif ~b1&&~b2
        lines2_attr(i,9) = -1;
    else
        lines2_attr(i,9) = 0;
    end   
end

%% save line segments according to type
typ = [1 1 2; 1 2 -1; 1 -1 2; 1 2 1; ...
       -1 1 2; -1 2 -1; -1 -1 2; -1 2 1; ...
       2 1 -1; 2 1 1; 2 -1 1; 2 -1 -1];
allLines = [lines1_attr; lines2_attr; lines3_attr];
typLine = repmat(struct('lines', [], 'lsLength', [], 'lspool', [], ...
                        'poolLength', []), [size(typ,1) 1]);
for i = 1:size(typ,1)
    valid = (allLines(:,7)==typ(i,1)) & (allLines(:,8)==typ(i,2)) & (allLines(:,9)==typ(i,3));
    typLine(i).lines = allLines(valid,:); 
    lsLength = typLine(i).lines(:,6)-typLine(i).lines(:,5);
    lsLength(lsLength<0) = lsLength(lsLength<0) + 1;
    sumLength = sum(lsLength);
    lspool = zeros( ceil(sumLength*100), 1);
    accumuLen = 0;
    poolLength = length(lspool);
    for j = 1:length(lsLength)
        startid = round(accumuLen/sumLength*poolLength) + 1;
        accumuLen = accumuLen + lsLength(j);
        endid = round(accumuLen/sumLength*poolLength);
        lspool(startid:endid) = j;
    end
    
    typLine(i).lsLength = lsLength;
    typLine(i).lspool = lspool;
    typLine(i).poolLength = poolLength;
end

%% >> BUILD UP RECIPE FROM LINE TYPE TO BOX
% +1: max coords; -1: min coords;
% uf [ 0 +1 +1]  df [ 0 +1 -1]  lf [-1 +1  0];
% ul [-1  0 +1]  dl [-1  0 -1]  rf [+1 +1  0];
% ub [ 0 -1 +1]  db [ 0 -1 -1]  rb [+1 -1  0];
% ur [+1  0 +1]  dr [+1  0 -1]  lb [-1 -1  0];

%% generate recipe: from line type to box by check rank of matrix
allFunc = [ 0 0 0 1 0 -1; ...
            1 0 0 0 0 -1; ...
            0 0 1 0 0 -1; ...
            0 1 0 0 0 -1; ...
            0 0 0 1 -1 0; ...
            1 0 0 0 -1 0; ...
            0 0 1 0 -1 0; ...
            0 1 0 0 -1 0; ...
            1 0 0 -1 0 0; ...
            0 1 0 -1 0 0; ...
            0 1 -1 0 0 0; ...
            1 0 -1 0 0 0];
% combos = combntns(1:size(allFunc,1), 5);
combos = nchoosek(1:size(allFunc,1), 5);
degen = false(size(combos,1),1);
for i = 1:size(combos,1)
    equset = allFunc(combos(i,:), :);
    if rank(equset)<5
        degen(i) = true;
    end
end
cube_recipe = combos(~degen,:);

%% define geometry of cuboid: correlation among corner, edge, and face
% 8 corners, each is the intersection of 3 edges from 3 vps respectively
edgeToPoint = [ 9  2  1; ...
               12  2  3; ...
               11  4  3; ...
               10  4  1; ...
                9  6  5; ...
               12  6  7; ...
               11  8  7; ...
               10  8  5];
% 8 corners, each is assigned with type [u(+1)/d(-1) f(+1)/b(-1) r(+1)/l(-1)]
pointDirect = [ 1  1 -1; ...
                1 -1 -1; ...
                1 -1  1; ...
                1  1  1; ...
               -1  1 -1; ...
               -1 -1 -1; ...
               -1 -1  1; ...
               -1  1  1];
% 12 edges, each connects two corners
PointToEdge = [ 1  4; ...
                2  1; ...
                2  3; ...
                3  4; ...
                5  8; ...
                6  5; ...
                6  7; ...
                7  8; ...
                5  1; ...
                8  4; ...
                7  3; ...
                6  2];
% 6 faces, each face has 4 corners
% ranked as u d f b r l; start from left bottom, clockwise
PointToSurface = [1 2 3 4; ...
                  6 5 8 7; ...
                  5 1 4 8; ...
                  7 3 2 6; ...
                  8 4 3 7; ...
                  6 2 1 5];
SurfaceNm = [1;1;2;2;3;3];

%% >> START TO GENERATE CUBOID HYPOTHESIS, VALIDATE WITH OMAP
maxSample = 200000;
validHyps = repmat(struct('extLine', [], 'hCorner', [], 'srcLine', [], ...
                          'omapScr', 0, 'vpLength', [], 'recipe', [], 'gcScr', -1, 'mgScr', -1, 'opScr', -1 ), maxSample, 1);

%% sample line segments, long line segments are more likely to be selected.
rcpNum = size(cube_recipe,1);
for samid = 1:maxSample
    if rem(samid,1000)==0
%         fprintf('%d/%d hypothesis validated, score: %f/%f\n', samid, maxSample, maxScore, numPoint);
        fprintf('%d/%d hypothesis sampled\r', samid, maxSample);
    end
    rcpid = cube_recipe(randsample(rcpNum, 1),:);
    
    if any([typLine(rcpid).poolLength]==0)
        validHyps(samid).srcLine = [];
        validHyps(samid).recipe = rcpid;  
        continue;
    end
    
    lineID1 = typLine(rcpid(1)).lspool(randsample(typLine(rcpid(1)).poolLength, 1));
    lineID2 = typLine(rcpid(2)).lspool(randsample(typLine(rcpid(2)).poolLength, 1));
    lineID3 = typLine(rcpid(3)).lspool(randsample(typLine(rcpid(3)).poolLength, 1));
    lineID4 = typLine(rcpid(4)).lspool(randsample(typLine(rcpid(4)).poolLength, 1));
    lineID5 = typLine(rcpid(5)).lspool(randsample(typLine(rcpid(5)).poolLength, 1));
    
    hlines = zeros(12,9);
    hlines(rcpid(1),:) = typLine(rcpid(1)).lines(lineID1,:);
    hlines(rcpid(2),:) = typLine(rcpid(2)).lines(lineID2,:);
    hlines(rcpid(3),:) = typLine(rcpid(3)).lines(lineID3,:);
    hlines(rcpid(4),:) = typLine(rcpid(4)).lines(lineID4,:);
    hlines(rcpid(5),:) = typLine(rcpid(5)).lines(lineID5,:);
    
    validHyps(samid).srcLine = hlines;
    validHyps(samid).recipe = rcpid;
end
fprintf('\n');

%% finish up unknown edges, compute corner, check if ls locates in line
for samid = 1:maxSample
    if rem(samid,1000)==0
%         fprintf('%d/%d hypothesis validated, score: %f/%f\n', samid, maxSample, maxScore, numPoint);
        fprintf('%d/%d hypothesis validated\r', samid, maxSample);
    end
    if isempty(validHyps(samid).srcLine)
        validHyps(samid).omapScr = -1;
        continue;
    end
    hlines = validHyps(samid).srcLine;
    hCorner = zeros(8,3);
    lineExist = false(12,1);
    lineExist(validHyps(samid).recipe) = true;
    rcpid = validHyps(samid).recipe;
    % >> compute all the other edges, corner
    bDegenerate = false;
    while any(~lineExist)
        hits = lineExist(edgeToPoint);
        actLine = sum(hits,2)==2;
        i = find(actLine,1)';
        srcLineID = edgeToPoint(i,hits(i,:));
        dstLineID = edgeToPoint(i,~hits(i,:));
        dstVpID = find(~hits(i,:));
        dstCor = cross( hlines(srcLineID(1),1:3), hlines(srcLineID(2),1:3), 2);
        if norm(dstCor)<0.1
            bDegenerate = true;
            break;
        end

        hCorner(i,:) = dstCor./norm(dstCor);
        % check if intersection locates at right position
        v = dot([dstCor;dstCor;dstCor],[horizontal;frontNback;leftNright],2);
        if sign(v(1))~=pointDirect(i,1) ...
         ||sign(v(2))~=pointDirect(i,2) ...
         ||sign(v(3))~=pointDirect(i,3)
            hCorner(i,:) = -hCorner(i,:);
        end

        dstNM = cross(hCorner(i,:), vp(dstVpID,:), 2);
        hlines(dstLineID,1:3) = dstNM./norm(dstNM);   
        lineExist(dstLineID) = true;
    end
    
    if bDegenerate
        validHyps(samid).omapScr = -1;
        continue;
    end
    
    I = find(actLine); i = I(end);
    srcLineID = edgeToPoint(i,hits(i,:));
    dstCor = cross( hlines(srcLineID(1),1:3), hlines(srcLineID(2),1:3), 2);
    hCorner(i,:) = dstCor./norm(dstCor);
    % check if intersection locates at right position
    v = dot([dstCor;dstCor;dstCor],[horizontal;frontNback;leftNright],2);
    if sign(v(1))~=pointDirect(i,1) ...
     ||sign(v(2))~=pointDirect(i,2) ...
     ||sign(v(3))~=pointDirect(i,3)
        hCorner(i,:) = -hCorner(i,:);
    end
    
    % >> compute extended line representation
    rcpidR = setdiff(1:12, rcpid);
    areaXY = abs(sum(hlines(rcpidR,1:3).*repmat([0 0 1], [length(rcpidR) 1]),2));
    areaYZ = abs(sum(hlines(rcpidR,1:3).*repmat([1 0 0], [length(rcpidR) 1]),2));
    areaZX = abs(sum(hlines(rcpidR,1:3).*repmat([0 1 0], [length(rcpidR) 1]),2));
    [~, planeIDs] = max([areaXY areaYZ areaZX], [], 2); % 1:XY 2:YZ 3:ZX
    hlines(rcpidR,4) = planeIDs;

    extLine = hlines(:,1:6);
    for i = 1:12
        endPoint = hCorner(PointToEdge(i,:),:);
        uv = xyz2uvN(endPoint,extLine(i,4));
        umin = (min(uv(:,1))+pi)/2/pi;    umax = (max(uv(:,1))+pi)/2/pi;
        if umax-umin>0.5
            extLine(i,5) = umax;
            extLine(i,6) = umin;
        else
            extLine(i,5) = umin;
            extLine(i,6) = umax;
        end
    end
    
    % check if line segments locates in line
    validated = true;
    for i = rcpid
        lineRange = extLine(i,5:6);
        lsegRange = hlines(i,5:6);
        b = insideRange(lsegRange, lineRange);
        if ~b(1) || ~b(2)
            validated = false;
            break;
        end
    end
    if ~validated
        validHyps(samid).omapScr = -1;
        continue;
    end
    validHyps(samid).srcLine = hlines(rcpid,1:6);
    validHyps(samid).extLine = extLine;
    validHyps(samid).hCorner = hCorner;
end
fprintf('\n');

%% compute the checking direction on omap                     
[candiSetXYZ, ~] = icosahedron2sphere(6);
[ohei, owid, ~] = size(omap); 
candiSetUV = uv2coords(xyz2uvN(candiSetXYZ), owid, ohei);
candiInd = sub2ind([ohei owid], candiSetUV(:,2), candiSetUV(:,1));
% convert gc to omap
gndGC = zeros(ohei, owid, 3);
gndGC(:,:,1) = max(gc(:,:,5),gc(:,:,6));
gndGC(:,:,2) = max(gc(:,:,1),gc(:,:,3));
gndGC(:,:,3) = max(gc(:,:,2),gc(:,:,4));
normGndGC = sum(gndGC(candiInd))+sum(gndGC(candiInd+1*ohei*owid))+sum(gndGC(candiInd+2*ohei*owid));
% omap
gndOmap = omap;
gndOmap(gndOmap>0) = 1;
normGndOmap = sum(gndOmap, 3);
gndOmap = gndOmap./(repmat(normGndOmap+0.0001, [1 1 3])); 
normGndOmap = sum(gndOmap(candiInd))+sum(gndOmap(candiInd+1*ohei*owid))+sum(gndOmap(candiInd+2*ohei*owid));
numPoint = size(candiSetXYZ,1);

% merge
gndMerge = zeros(ohei, owid, 3);
gndMerge(1:ohei/2,:,:) = gndOmap(1:ohei/2,:,:);
gndMerge(ohei/2+1:ohei, :, :) = gndGC(ohei/2+1:ohei, :, :);
% normMegMap = sum(gndOmap, 3);
% gndOmap = gndOmap./(repmat(normMegMap+0.0001, [1 1 3])); 
normMegMap = sum(gndMerge(candiInd))+sum(gndMerge(candiInd+1*ohei*owid))+sum(gndMerge(candiInd+2*ohei*owid));

% optimal merge
gndOptim = zeros(ohei, owid, 3);
gndOptim(1:685,:,:) = gndOmap(1:685,:,:);
gndOptim(686:ohei, :, :) = gndGC(686:ohei, :, :);
% normMegMap = sum(gndOmap, 3);
% gndOmap = gndOmap./(repmat(normMegMap+0.0001, [1 1 3])); 
normOptMap = sum(gndOptim(candiInd))+sum(gndOptim(candiInd+1*ohei*owid))+sum(gndOptim(candiInd+2*ohei*owid));

%% check consistency of orientation on uniformly sampled point on sphere  
for samid = 1:maxSample
    hypVp = false(numPoint,3);
    if validHyps(samid).omapScr>-0.5
        hCorner = validHyps(samid).hCorner;
        for i = 1:6
%             suf1 = cross(hCorner(PointToSurface(i,2),:), hCorner(PointToSurface(i,1),:), 2);
%             suf2 = cross(hCorner(PointToSurface(i,3),:), hCorner(PointToSurface(i,2),:), 2);
%             suf3 = cross(hCorner(PointToSurface(i,4),:), hCorner(PointToSurface(i,3),:), 2);
%             suf4 = cross(hCorner(PointToSurface(i,1),:), hCorner(PointToSurface(i,4),:), 2);
%             vad1 = dot(candiSetXYZ, repmat(suf1, [numPoint 1]), 2);
%             vad2 = dot(candiSetXYZ, repmat(suf2, [numPoint 1]), 2);
%             vad3 = dot(candiSetXYZ, repmat(suf3, [numPoint 1]), 2);
%             vad4 = dot(candiSetXYZ, repmat(suf4, [numPoint 1]), 2);
%             vadIN = vad1<0 & vad2<0 & vad3<0 & vad4<0;
            vadIN = insideCone( hCorner(PointToSurface(i,[4 3 2 1]),:), candiSetXYZ, 0);
            hypVp( vadIN, SurfaceNm(i)) = true;
        end
%         response = selVp & hypVp;
%         validHyps(samid).omapScr = sum( sum(response,2)>0 )/size(response,1);
        response = sum(gndOmap(candiInd(hypVp(:,1))+0*ohei*owid)) ...
                 + sum(gndOmap(candiInd(hypVp(:,2))+1*ohei*owid)) ...
                 + sum(gndOmap(candiInd(hypVp(:,3))+2*ohei*owid));
        validHyps(samid).omapScr = response/normGndOmap;
        response = sum(gndGC(candiInd(hypVp(:,1))+0*ohei*owid)) ...
                 + sum(gndGC(candiInd(hypVp(:,2))+1*ohei*owid)) ...
                 + sum(gndGC(candiInd(hypVp(:,3))+2*ohei*owid));
        validHyps(samid).gcScr = response/normGndGC;
        response = sum(gndMerge(candiInd(hypVp(:,1))+0*ohei*owid)) ...
                 + sum(gndMerge(candiInd(hypVp(:,2))+1*ohei*owid)) ...
                 + sum(gndMerge(candiInd(hypVp(:,3))+2*ohei*owid));
        validHyps(samid).mgScr = response/normMegMap;
        response = sum(gndOptim(candiInd(hypVp(:,1))+0*ohei*owid)) ...
                 + sum(gndOptim(candiInd(hypVp(:,2))+1*ohei*owid)) ...
                 + sum(gndOptim(candiInd(hypVp(:,3))+2*ohei*owid));
        validHyps(samid).opScr = response/normOptMap;
    end
end

%% save only validate hypothesis, with high score
omapScr = [validHyps.omapScr];
% [B, IX] = sort(omapScr, 'descend');
% selIX = IX(B>-0.5);
% hyps = validHyps(selIX);

hyps = validHyps(omapScr>-0.5);
end

function b = insideRange(pt, range)
range = range + [-0.02 +0.02];
if range(2)>range(1)
    b = pt>=range(1) & pt<=range(2);
else
    b1 = pt>=range(1) & pt<=1;
    b2 = pt>=0 & pt<=range(2);
    b = b1 | b2;
end

end
