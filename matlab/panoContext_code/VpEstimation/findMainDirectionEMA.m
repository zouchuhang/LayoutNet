function [ mainDirect, score, angle] = findMainDirectionEMA( lines )
%FINDMAINDIRECTION compute vp from set of lines
%   Detailed explanation goes here
fprintf('Computing vanishing point:\n');

% arcList = [];
% for i = 1:length(edge)
%     panoLst = edge(i).panoLst;
%     if size(panoLst,1) == 0
%         continue;
%     end
%     arcList = [arcList; panoLst];
% end

%% initial guess
segNormal = lines(:,1:3);
segLength = lines(:,7);
segScores = ones(size(lines,1),1);%lines(:,8);

shortSegValid = segLength < 5*pi/180;
segNormal = segNormal(~shortSegValid,:);
segLength = segLength(~shortSegValid);
segScores = segScores(~shortSegValid);

numLinesg = size(segNormal,1);
[candiSet, tri] = icosahedron2sphere(3);
ang = acos(dot(candiSet(tri(1,1),:), candiSet(tri(1,2),:), 2)) / pi * 180;
binRadius = ang/2;
[ initXYZ, score, angle] = sphereHoughVote( segNormal, segLength, segScores, 2*binRadius, 2, candiSet );

if isempty(initXYZ)
    fprintf('Initial Failed\n');
    mainDirect = [];
    return;
end

fprintf('Initial Computation: %d candidates, %d line segments\n', size(candiSet,1), numLinesg);
fprintf('direction 1: %f %f %f\ndirection 2: %f %f %f\ndirection 3: %f %f %f\n', ...
        initXYZ(1,1), initXYZ(1,2), initXYZ(1,3), ...
        initXYZ(2,1), initXYZ(2,2), initXYZ(2,3), ...
        initXYZ(3,1), initXYZ(3,2), initXYZ(3,3));
%% iterative refine
iter_max = 3;
[candiSet, tri] = icosahedron2sphere(5);
numCandi = size(candiSet,1);
angD = acos(dot(candiSet(tri(1,1),:), candiSet(tri(1,2),:), 2)) / pi * 180;
binRadiusD = angD/2;
curXYZ = initXYZ;
tol = linspace(4*binRadius, 4*binRadiusD, iter_max); % shrink down #ls and #candi
for iter = 1:iter_max
    dot1 = abs(dot( segNormal, repmat(curXYZ(1,:), [numLinesg 1]), 2));
    dot2 = abs(dot( segNormal, repmat(curXYZ(2,:), [numLinesg 1]), 2));
    dot3 = abs(dot( segNormal, repmat(curXYZ(3,:), [numLinesg 1]), 2));
    valid1 = dot1<cos((90-tol(iter))*pi/180);
    valid2 = dot2<cos((90-tol(iter))*pi/180);
    valid3 = dot3<cos((90-tol(iter))*pi/180);
    valid = valid1 | valid2 | valid3;
    
    if(sum(valid)==0)
        fprintf('ZERO line segment for voting\n');
        break;
    end
    
    subSegNormal = segNormal(valid,:);
    subSegLength = segLength(valid);
    subSegScores = segScores(valid);
    
    dot1 = abs(dot( candiSet, repmat(curXYZ(1,:), [numCandi 1]), 2));
    dot2 = abs(dot( candiSet, repmat(curXYZ(2,:), [numCandi 1]), 2));
    dot3 = abs(dot( candiSet, repmat(curXYZ(3,:), [numCandi 1]), 2));
    valid1 = dot1>cos(tol(iter)*pi/180);
    valid2 = dot2>cos(tol(iter)*pi/180);
    valid3 = dot3>cos(tol(iter)*pi/180);
    valid = valid1 | valid2 | valid3;
    
    if(sum(valid)==0)
        fprintf('ZERO candidate for voting\n');
        break;
    end
       
    subCandiSet = candiSet(valid,:);
    
    [ tcurXYZ ] = sphereHoughVote( subSegNormal, subSegLength, subSegScores, 2*binRadiusD, 2, subCandiSet );
    
    if(isempty(tcurXYZ))
        fprintf('NO answer found!\n');
        break;
    end
    curXYZ = tcurXYZ;

    fprintf('%d-th iteration: %d candidates, %d line segments\n', iter, size(subCandiSet,1), length(subSegScores));

end
fprintf('direction 1: %f %f %f\ndirection 2: %f %f %f\ndirection 3: %f %f %f\n', ...
        curXYZ(1,1), curXYZ(1,2), curXYZ(1,3), ...
        curXYZ(2,1), curXYZ(2,2), curXYZ(2,3), ...
        curXYZ(3,1), curXYZ(3,2), curXYZ(3,3));
mainDirect = curXYZ;

mainDirect(1,:) = mainDirect(1,:).*sign(mainDirect(1,3));
mainDirect(2,:) = mainDirect(2,:).*sign(mainDirect(2,3));
mainDirect(3,:) = mainDirect(3,:).*sign(mainDirect(3,3));

uv = xyz2uvN(mainDirect);
[~,I1] = max(uv(:,2));
J = setdiff(1:3, I1);
[~,I2] = min(abs(sin(uv(J,1))));
I2 = J(I2);
I3 = setdiff(1:3, [I1 I2]);
mainDirect = [mainDirect(I1,:); mainDirect(I2,:); mainDirect(I3,:)];

mainDirect(1,:) = mainDirect(1,:)*sign(mainDirect(1,3));
mainDirect(2,:) = mainDirect(2,:)*sign(mainDirect(2,2));
mainDirect(3,:) = mainDirect(3,:)*sign(mainDirect(3,1));

mainDirect = [mainDirect; -mainDirect];


% score = 0;

end

