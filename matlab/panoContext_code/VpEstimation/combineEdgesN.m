function [ lines, ori_lines ] = combineEdgesN(  edge )
%COMBINEEDGES Combine some small line segments, should be very conservative
%   lines: combined line segments
%   ori_lines: original line segments
%   line format: [nx ny nz projectPlaneID umin umax LSfov score]
arcList = [];
for i = 1:length(edge)
    panoLst = edge(i).panoLst;
    if size(panoLst,1) == 0
        continue;
    end
    arcList = [arcList; panoLst];
end

%% ori lines
numLine = size(arcList,1);
ori_lines = zeros(numLine,8);
areaXY = abs(sum(arcList(:,1:3).*repmat([0 0 1], [size(arcList,1) 1]),2));
areaYZ = abs(sum(arcList(:,1:3).*repmat([1 0 0], [size(arcList,1) 1]),2));
areaZX = abs(sum(arcList(:,1:3).*repmat([0 1 0], [size(arcList,1) 1]),2));
[~, planeIDs] = max([areaXY areaYZ areaZX], [], 2); % 1:XY 2:YZ 3:ZX

for i = 1:numLine
    ori_lines(i,1:3) = arcList(i,1:3);
    ori_lines(i,4) = planeIDs(i);
    coord1 = arcList(i,4:6);
    coord2 = arcList(i,7:9);
    uv = xyz2uvN([coord1; coord2], planeIDs(i));
    umax = max(uv(:,1))+pi;
    umin = min(uv(:,1))+pi;
    if umax-umin>pi
        ori_lines(i,5:6) = [umax umin]/2/pi;
%         ori_lines(i,7) = umin + 1 - umax;
    else
        ori_lines(i,5:6) = [umin umax]/2/pi;
%         ori_lines(i,7) = umax - umin;
    end
    ori_lines(i,7) = acos(dot(coord1, coord2, 2)./(norm(coord1)*norm(coord2)));
    ori_lines(i,8) = arcList(i,10);
end
% valid = ori_lines(:,3)<0;
% ori_lines(valid,1:3) = -ori_lines(valid,1:3);


%% additive combination
lines = ori_lines;
% panoEdge = paintParameterLine( lines, 1024, 512);
% figure; imshow(panoEdge);
for iter = 1:3
    numLine = size(lines,1);
    valid_line = true(numLine,1);
    for i = 1:numLine
%         fprintf('%d/%d\n', i, numLine);
        if ~valid_line(i)
            continue;
        end
        dotProd = dot(lines(:,1:3), repmat(lines(i,1:3), [numLine 1]), 2);
        valid_curr =  (abs(dotProd) > cos(1*pi/180)) & valid_line;
        valid_curr(i) = false;
        valid_ang = find(valid_curr);
        for j = valid_ang'
            range1 = lines(i,5:6);
            range2 = lines(j,5:6);
            valid_rag = intersection(range1, range2);
            if ~valid_rag
                continue;
            end

            % combine   
            [~,I] = max(abs(lines(i,1:3)));
            if lines(i,I)*lines(j,I)>0
                nc = lines(i,1:3)*lines(i, 7) + lines(j,1:3)*lines(j, 7);
            else
                nc = lines(i,1:3)*lines(i, 7) - lines(j,1:3)*lines(j, 7);
            end
            nc = nc / norm(nc);

            if insideRange(range1(1), range2)
                nrmin = range2(1);
            else
                nrmin = range1(1);
            end
            if insideRange(range1(2), range2)
                nrmax = range2(2);
            else
                nrmax = range1(2);
            end

            u = [nrmin;nrmax]*2*pi - pi;
            v = computeUVN( nc, u, lines(i,4));
            xyz = uv2xyzN([u v], lines(i,4));
            len = acos(dot(xyz(1,:), xyz(2,:), 2));
            scr = (lines(i,7)*lines(i,8) + lines(j,7)*lines(j,8))/(lines(i,7)+lines(j,7));
            
            newLine = [nc lines(i,4) nrmin nrmax len scr];
            lines(i,:) = newLine;
            valid_line(j) = false;
        end
    end
    lines(~valid_line,:) = []; 
    fprintf('iter: %d, before: %d, after: %d\n', iter, length(valid_line), sum(valid_line));
%     panoEdge = paintParameterLine( lines, 1024, 512);
%     figure; imshow(panoEdge);
end
%% previous method, bin voting
% %% build up voting space
% numDivision = 15;%round(max(width, height)/4);
% 
% normal = arcList(:,1:3);
% valid = normal(:,3)<0;
% normal(valid,:) = -normal(valid,:);
% uv = xyz2uvN(normal, 1);
% thetas = uv(:,1);   phis = uv(:,2);
% 
% uBinSize = 2*pi/(4*numDivision);
% vBinSize = pi/2/(numDivision);
% m = min(floor( (thetas-(-pi))/uBinSize) + 1, numDivision*4);
% n = min(floor( phis/vBinSize) + 1, numDivision);
% normalInd = sub2ind([4*numDivision numDivision], m, n);
% 
% uniqueNormal = unique(normalInd);
% % decide voting dimension
% [m,n] = ind2sub([4*numDivision numDivision], uniqueNormal);
% u = -pi + (m-1)*uBinSize + uBinSize/2;
% v = 0   + (n-1)*vBinSize + vBinSize/2;
% uniNormal = [cos(v).*sin(u) cos(v).*cos(u) sin(v)];
% areaXY = abs(sum(uniNormal.*repmat([0 0 1], [size(uniNormal,1) 1]),2));
% areaYZ = abs(sum(uniNormal.*repmat([1 0 0], [size(uniNormal,1) 1]),2));
% areaZX = abs(sum(uniNormal.*repmat([0 1 0], [size(uniNormal,1) 1]),2));
% [~, planeIDs] = max([areaXY areaYZ areaZX], [], 2); % 1:XY 2:YZ 3:ZX
% 
% subVoteBinNum = 1024;
% uvBin = false(length(uniqueNormal), subVoteBinNum);
% uBinSize = 2*pi/subVoteBinNum;
% for i = 1:length(uniqueNormal)
%     normIDs = find(normalInd==uniqueNormal(i));
%     norms = normal(normIDs,:);
%     rectNorm = sum(norms,1);
%     uniNormal(i,:) = rectNorm./norm(rectNorm);
%     dimIndc = planeIDs(i);
%     for j = normIDs'
%         coord1 = arcList(j,4:6);
%         coord2 = arcList(j,7:9);
%         xx = linspace(coord1(1),coord2(1),1000);
%         yy = linspace(coord1(2),coord2(2),1000);
%         zz = linspace(coord1(3),coord2(3),1000);
%         uv = xyz2uvN([xx' yy' zz'], dimIndc);  
% 
%         m = min(floor( (uv(:,1)-(-pi))/uBinSize) + 1, subVoteBinNum);
%         uvBin(i,m) = true;
%     end
% end
% 
% %% extract line segments
% % numLines = 0;
% lines = [];
% for i = 1:length(uniqueNormal)
% %     fprintf('%d\n',i);
%     bins = int32(uvBin(i,:));
%     changePt(2:length(bins)) = bins(2:end) - bins(1:end-1);
%     changePt(1) = bins(1) - bins(end);
%     startPt = find(changePt==1);
%     endPt   = find(changePt==-1) - 1;
%     endPt = rem(endPt + subVoteBinNum+1, subVoteBinNum+1);
%     
%     if endPt(1)>=startPt(1)
%         mtEndPt = endPt;
%     else
%         mtEndPt = [endPt(2:end) endPt(1)]; 
%     end
%     
%     lines(end+1:end+length(startPt),:) = ...
%         [repmat([uniNormal(i,:) planeIDs(i)],[length(startPt) 1]) ...
%         (startPt'-1)/subVoteBinNum mtEndPt'/subVoteBinNum];
% end
end

function b = intersection(range1, range2)
    if range1(2)<range1(1)
        range11 = [range1(1) 1];
        range12 = [0 range1(2)];
    else
        range11 = range1;
        range12 = [0 0];
    end
    if range2(2)<range2(1)
        range21 = [range2(1) 1];
        range22 = [0 range2(2)];
    else
        range21 = range2;
        range22 = [0 0];
    end
    b = max(range11(1),range21(1))<min(range11(2),range21(2));
    if b
        return;
    end
    b2 = max(range12(1),range22(1))<min(range12(2),range22(2));
    b = b||b2;
end

function b = insideRange(pt, range)
if range(2)>range(1)
    b = pt>=range(1) && pt<=range(2);
else
    b1 = pt>=range(1) && pt<=1;
    b2 = pt>=0 && pt<=range(2);
    b = b1 || b2;
end

end
