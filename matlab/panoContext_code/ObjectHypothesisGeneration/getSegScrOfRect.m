function [ segment_rect_score ] = getSegScrOfRect( rotImg, rectangles )
%GETSEGSCROFRECT Refine rectangle detected by detector with segmentation
%   Each rectangle get a score as intersection/union with segmentation

%% compute multiple layer of segmentation
rotImg  = im2double(rotImg);
thresholds = [200 500 800 1200 2000 3000];
for i = 1:length(thresholds)
    all_seg(i).labelMap = gbPanoSegment( im2uint8(rotImg), 0.5, thresholds(i), 50);
    all_seg(i).thresh = thresholds(i);
end

%% checkseg
% load('./rectangleDetector/segmentation/uniformvector_lvl6.mat');
[coor, tri] = getUniformVector(6);
[sH, sW, ~] = size(all_seg(1).labelMap);
coor2D = uv2coords(xyz2uvN(coor), sW, sH);
coor2D_indi = sub2ind([sH sW], coor2D(:,2), coor2D(:,1));

% vp = [-1  0  0; ...
%        1  0  0; ...
%        0 -1  0; ...
%        0  1  0; ...
%        0  0 -1; ...
%        0  0  1];

for rid = 1:6
    xyzBox = rectangles(rid).xyzBox;
    vector_num = size(coor, 1);
    indicator = false(vector_num, size(xyzBox,1));
    for bid = 1:size(xyzBox,1)
        p = reshape( xyzBox(bid,:), 3, 4)';
        p = p./repmat(sqrt(sum(p.^2,2)),1,3);
        vp = sum(p,1); vp = vp./norm(vp); %vp(3) = 0; 
        [out2D, valid, ~] = projectPoint2SeparateView( p, vp, pi/3, 100 );
        if any(~valid)
            vp = sum(p,1); vp(1:2) = 0; vp = vp./norm(vp); 
            [out2D, valid, ~] = projectPoint2SeparateView( p, vp, pi/3, 100 );
            if any(~valid)
                vp = sum(p,1); vp(3) = 0; vp = vp./norm(vp); 
                [out2D, valid, ~] = projectPoint2SeparateView( p, vp, pi/3, 100 );
            end
        end
        try
            K = convhull( out2D(:,1), out2D(:,2));
            [ inside, ~, ~ ] = insideCone( p(K(end:-1:2),:), coor, 0 );
    %         if ~any(inside)
    %             fprintf('%d:%d\n', bid, i);
    %         end
    %         assert(any(inside));
            indicator(:,bid) = inside;
        catch
            continue;
        end
    end
    
    seg_score = zeros( length(all_seg), size(xyzBox,1));
    for sid = 1:length(all_seg)
        label = all_seg(sid).labelMap(coor2D_indi);
        MINLAB = min(label);
        MAXLAB = max(label);
        N = hist(label, MINLAB:MAXLAB);
        for j = 1:size(xyzBox,1)
            lab = label(indicator(:,j));
            if isempty(lab)
                continue;
            end
            M = hist(lab, MINLAB:MAXLAB);
            I = M./sum(M)>0.05 & M./N>0.05;
%             I = M>0.05*sum(M);
            if sum(I)==0
                seg_score(sid,j) = 0;
            else
                seg_score(sid,j) = sum(M(I))./max(sum(N(I)),sum(M));
            end
        end
    end
    max_seg_score = max(seg_score, [], 1);
    segment_rect_score{rid} = max_seg_score;
    
end

end

