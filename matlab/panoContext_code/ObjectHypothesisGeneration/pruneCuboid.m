function [ AllCuboid ] = pruneCuboid( rotImg, AllCuboid )
%PRUNECUBOID Prune cuboid by checking gradient on edges
%   Edges of cuboids should present relatively strong gradient

rotImg = imresize(rotImg, 0.25);
[H,W,~] = size(rotImg);
% img = zeros(H,W);

[FX, FY] = gradient(double(rgb2gray(rotImg)));
F = sqrt(FX.^2 + FY.^2);
[FD, ~, ~] = dt2(F, 0.1, 0, 0.1, 0 );

valid = false(AllCuboid.count, 1);

for cid = 1:AllCuboid.count
    view = AllCuboid.views(cid,:);
    xyzBox = AllCuboid.xyzBox(cid,:);
    numRectangle = sum(view~=0);
    AllLines = zeros(0,6);
    for rid = 1:numRectangle
        rect = reshape(xyzBox((rid*12-11):(rid*12)), 3, 4)';
        lines = lineFromTwoPoint(rect([1 2 3 4], :), rect([2 3 4 1], :));
        AllLines = [AllLines; lines];
    end
    panoReport = paintParameterLine( AllLines, W, H);
    G = FD(panoReport(:)>0);
    SG = sort(G, 'descend');
    N = round(length(G)/2);

%         valid(i) = sum(SG(1:N))/N;
    if sum(SG(1:N))/N >20
        valid(cid) = true;
    end
    
end

AllCuboid.views = AllCuboid.views(valid,:);
AllCuboid.xyzBox = AllCuboid.xyzBox(valid,:);
AllCuboid.score = AllCuboid.score(valid,:);
AllCuboid.count = sum(valid);


end

