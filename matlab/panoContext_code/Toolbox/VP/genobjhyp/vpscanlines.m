function scanlineendpts = vpscanlines(vp, scandir, imgsize, OMAP_FACTOR)
currvp = OMAP_FACTOR * vp{scandir};

% if vp inside image
if currvp(1)>=1 && currvp(1)<=imgsize(2) && ...
   currvp(2)>=1 && currvp(2)<=imgsize(1)
    pivotpts = [ [(1:imgsize(2))' repmat(1,imgsize(2),1)]; ... %top
                 [(1:imgsize(2))' repmat(imgsize(1),imgsize(2),1)]; ...%bottom
                 [repmat(imgsize(2),imgsize(1),1) (1:imgsize(1))']; ... %right
                 [repmat(        1, imgsize(1),1) (1:imgsize(1))'] ]; %left

    scanlineendpts = zeros(size(pivotpts,1), 4);
    scanlineendpts(:,1:2) = pivotpts;
    scanlineendpts(:,3:4) = repmat(currvp, size(pivotpts,1), 1);

% if outside image
else
    pivotpts = [];
    if currvp(1) < 1
        pivotpts = [repmat(         1,imgsize(1),1) (1:imgsize(1))'];
    elseif currvp(1) > imgsize(2)
        pivotpts = [repmat(imgsize(2),imgsize(1),1) (1:imgsize(1))'];
    end

    if currvp(2) < 1
        pivotpts = [pivotpts; ...
                    (1:imgsize(2))' repmat(1,imgsize(2),1)];
    elseif currvp(2) > imgsize(1)
        pivotpts = [pivotpts; ...
                    (1:imgsize(2))' repmat(imgsize(1),imgsize(2),1)];
    end

    scanlineendpts = zeros(size(pivotpts,1), 4);
    for i = 1:size(pivotpts,1)
        [p1 p2] = extline(pivotpts(i,:), currvp, imgsize(2), imgsize(1));
        
        % TODO: change here. change from distance to camera to up/down, left/right
        % DONE: I think the above TODO has been done.. 05/22/2010 --dcl
        if scandir==1
            if p1(2) > p2(2)
                temp = p1;
                p1 = p2;
                p2 = temp;
            end
        elseif scandir==2
            if p1(1) > p2(1)
                temp = p1;
                p1 = p2;
                p2 = temp;
            end
        elseif scandir==3
            if norm(p1-currvp) < norm(p2-currvp)
                temp = p1;
                p1 = p2;
                p2 = temp;
            end % p2 should be closer to vp than p1
        end
        
        scanlineendpts(i,1:2) = p1;
        scanlineendpts(i,3:4) = p2;
    end
end

% drop scanlines that are too short
scanlinelength = sqrt( ...
    (scanlineendpts(:,1)-scanlineendpts(:,3)).^2 + ...
    (scanlineendpts(:,2)-scanlineendpts(:,4)).^2 );
scanlineendpts = scanlineendpts(scanlinelength>5,:);


% if scandir==1
%     if currvp(2) < 0 % above image
%         pivotpts = [(1:imgsize(2))' repmat(1,imgsize(2),1)];
%     else % below image
%         pivotpts = [(1:imgsize(2))' repmat(imgsize(1),imgsize(2),1)];
%     end
%     if size(pivotpts,1) > 10
%         pivotpts = pivotpts(6:end-5,:); %remove those too close to corner
%     end
% elseif scandir==2
%     if currvp(1) < 0 % left of image
%         pivotpts = [repmat(         1,imgsize(1),1) (1:imgsize(1))'];
%     else % right of image
%         pivotpts = [repmat(imgsize(2),imgsize(1),1) (1:imgsize(1))'];
%     end
%     if size(pivotpts,1) > 10
%         pivotpts = pivotpts(6:end-5,:); %remove those too close to corner
%     end
% elseif scandir==3
%     pivotpts = [ [(1:imgsize(2))' repmat(1,imgsize(2),1)]; ... %top
%                  [(1:imgsize(2))' repmat(imgsize(1),imgsize(2),1)]; ...%bottom
%                  [repmat(imgsize(2),imgsize(1),1) (1:imgsize(1))']; ... %right
%                  [repmat(        1, imgsize(1),1) (1:imgsize(1))'] ]; %left
% end

