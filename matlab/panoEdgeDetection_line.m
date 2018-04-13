function [ olines, vp, views, edges, panoEdge, score, angle] ...
            = panoEdgeDetection( img, viewSize, qError )
%PANOEDGEDETECTION line detection on panorama
%   INPUT:
%   img: image waiting for detection, double type, range 0~1
%   viewSize: image size of croped views
%   qError: set smaller if more line segment wanted
%   OUTPUT:
%   oLines: detected line segments
%   vp: vanishing point
%   views: separate views of panorama
%   edges: original detection of line segments in separate views
%   panoEdge: image for visualize line segments

cutSize = viewSize;
fov = pi/3;
xh = -pi:(pi/6):(5/6*pi);
yh = zeros(1, length(xh));
xp = [-3/3 -2/3 -1/3 +0/3 +1/3 +2/3 -3/3 -2/3 -1/3 +0/3 +1/3 +2/3] * pi;
yp = [ 1/4  1/4  1/4  1/4  1/4  1/4 -1/4 -1/4 -1/4 -1/4 -1/4 -1/4] * pi;
x = [xh xp 0     0];
y = [yh yp +pi/2 -pi/2];

%% line segment detection
% [eout,thresh] = edge_m(rgb2gray(img), 'canny');
[sepScene] = separatePano( img, fov, x, y, cutSize);
% for i = 1:length(sepScene)
%     sepScene(i).vx = rem(sepScene(i).vx + pi/6 + pi, 2*pi) - pi;
% end
% for i = 1:length(sepScene)
%     sepScene(i).img = tsmooth(sepScene(i).img, 0.005, 3, 0.01);
% end

numScene = length(sepScene);
edge(numScene) = struct('img',[],'edgeLst',[],'vx',[],'vy',[],'fov',[]);
for i = 1:numScene
    cmdline = sprintf('-q %f ', qError);
    [ edgeMap, edgeList ] = lsdWrap( sepScene(i).img, cmdline);
%     [edgeList, edgeMap] = getLargeConnectedEdges(sepScene(i).img, thresh);
    edge(i).img = edgeMap;
    edge(i).edgeLst = edgeList;
    edge(i).fov = sepScene(i).fov;
    edge(i).vx = sepScene(i).vx;
    edge(i).vy = sepScene(i).vy;
    edge(i).panoLst = edgeFromImg2Pano( edge(i) );
end
% [~,lines] = combineEdgesN( edge);
[lines,olines] = combineEdgesN( edge);
% panoEdge = paintParameterLine( lines, 1024, 512);
% figure; imshow(panoEdge);

%% compute vanishing point and refine line segments
% lines = olines;
% [mainDirect, ~] = findMainDirectionEM( edge );
clines = lines;
for iter = 1:3
    fprintf('*************%d-th iteration:****************\n', iter);
    [mainDirect, score, angle] = findMainDirectionEMA( clines );
    
%     [ type, typeCost ] = assignVanishingType( lines, mainDirect(1:3,:), 0.1, 10 );  
%     clines = lines(type<=3,:);
    
    
    [ type, typeCost ] = assignVanishingType( lines, mainDirect(1:3,:), 0.1, 10 ); 
%     lines1 = getVanishingLine(clines, mainDirect(1,:), 6);
%     lines2 = getVanishingLine(clines, mainDirect(2,:), 6);
%     lines3 = getVanishingLine(clines, mainDirect(3,:), 6);
    lines1 = lines(type==1,:);
    lines2 = lines(type==2,:);
    lines3 = lines(type==3,:);

%     panoEdge1 = paintParameterLine( lines1, 1024, 512);
%     panoEdge2 = paintParameterLine( lines2, 1024, 512);
%     panoEdge3 = paintParameterLine( lines3, 1024, 512);
%     panoEdge = cat(3, panoEdge1, panoEdge2, panoEdge3);
%     figure; subplot(2,1,1); imshow(panoEdge);
    
    lines1rB = refitLineSegmentB(lines1, mainDirect(1,:), 0);
    lines2rB = refitLineSegmentB(lines2, mainDirect(2,:), 0);
    lines3rB = refitLineSegmentB(lines3, mainDirect(3,:), 0);
    
%     panoEdge1r = paintParameterLine( lines1rB, 1024, 512);
%     panoEdge2r = paintParameterLine( lines2rB, 1024, 512);
%     panoEdge3r = paintParameterLine( lines3rB, 1024, 512);
%     panoEdger = cat(3, panoEdge1r, panoEdge2r, panoEdge3r);
%     subplot(2,1,2); imshow(panoEdger);
    
    clines = [lines1rB;lines2rB;lines3rB];
end

% [ type, typeCost ] = assignVanishingType( lines, mainDirect(1:3,:), 0.1, 10 );
% lines1rB = lines(type==1,:);
% lines2rB = lines(type==2,:);
% lines3rB = lines(type==3,:);
% clines = [lines1rB;lines2rB;lines3rB];

%imgres = imresize(img, [512 1024]);
imgres = zeros(size(img));
panoEdge1r = paintParameterLine( lines1rB, 1024, 512, imgres);
panoEdge2r = paintParameterLine( lines2rB, 1024, 512, imgres);
panoEdge3r = paintParameterLine( lines3rB, 1024, 512, imgres);
panoEdger = cat(3, panoEdge1r, panoEdge2r, panoEdge3r);
% figure; imshow(panoEdger);
% viewVanishingPoint(panoEdger, mainDirect);

%% output
olines = clines;
vp = mainDirect;

views = sepScene;
edges = edge;
panoEdge = panoEdger;

%% generate Hypothesis and validate

% [ omap, panoOmap ] = computePanoOmap( sepScene, edge, mainDirect );
% hyps = generateHyps(clines, mainDirect, 6, panoOmap);
% viewCorner(img, hyps);



% lines1rA = refitLineSegmentA(lines1, mainDirect(1,:));
% lines2rA = refitLineSegmentA(lines2, mainDirect(2,:));
% lines3rA = refitLineSegmentA(lines3, mainDirect(3,:));
% panoEdge1r = paintParameterLine( lines1rA, 1024, 512);
% panoEdge2r = paintParameterLine( lines2rA, 1024, 512);
% panoEdge3r = paintParameterLine( lines3rA, 1024, 512);
% panoEdger = cat(3, panoEdge1r, panoEdge2r, panoEdge3r);
% % imwrite(panoEdger, 'panoEdger.png');
% figure; imshow(panoEdger);
% 
% lines1rB = refitLineSegmentB(lines1, mainDirect(1,:));
% lines2rB = refitLineSegmentB(lines2, mainDirect(2,:));
% lines3rB = refitLineSegmentB(lines3, mainDirect(3,:));
% panoEdge1r = paintParameterLine( lines1rB, 1024, 512);
% panoEdge2r = paintParameterLine( lines2rB, 1024, 512);
% panoEdge3r = paintParameterLine( lines3rB, 1024, 512);
% panoEdger = cat(3, panoEdge1r, panoEdge2r, panoEdge3r);
% % imwrite(panoEdger, 'panoEdgerA.png');
% figure; imshow(panoEdger);

% [mainDirect, ~] = findMainDirectionEM( edge );
% lines1 = getVanishingLine(lines, mainDirect(1,:), 6);
% lines1 = lines1([1110 1133 1262],:);
% lines1r = refitLineSegment(lines1, mainDirect(1,:));
% panoEdge1 = paintParameterLine( lines1, 1024, 512);
% panoEdge1r = paintParameterLine( lines1r, 1024, 512);
% figure; imshow(panoEdge1);
% figure; imshow(panoEdge1r);
%% align lines to vanishing point


% [mainDirect, score] = findMainDirection( edge );




%% texture removal
% panoEdgeSm = panoEdge;
% sepSceneSm = sepScene;
% edgeSm = edge;
% sepSceneSm = sepScene;
% for i = 1:length(sepScene)
%     sepSceneSm(i).img = tsmooth(sepScene(i).img, 0.005, 3, 0.01);
% end
% edgeSm(numScene) = struct('img',[],'edgeLst',[],'vx',[],'vy',[],'fov',[]);
% for i = 1:numScene
%     cmdline = sprintf('-q %f', qError);
%     [ edgeMap, edgeList ] = lsdWrap( sepSceneSm(i).img, cmdline);
%     edgeSm(i).img = edgeMap;
%     edgeSm(i).edgeLst = edgeList;
%     edgeSm(i).fov = sepScene(i).fov;
%     edgeSm(i).vx = sepScene(i).vx;
%     edgeSm(i).vy = sepScene(i).vy;
%     edgeSm(i).panoLst = edgeFromImg2Pano( edgeSm(i) );
% end
% panoEdgeSm = combineViews( edgeSm, 1024, 512 );
% panoEdgeSm(panoEdgeSm>0) = 1;

% [lines, nBins, panoEdgeC] = combineEdgesN(  edge(1), 1024, 512  );

% figure(1);
% subplot(2,1,1); imshow(panoEdge);
% subplot(2,1,2); imshow(panoEdgeC);
end

