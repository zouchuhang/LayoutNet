function [v2 vpvar p hpos]=getVP(All_lines,imsize,DO_DISPLAY) 

    addpath('E:\Research\FunctionalLabeling\GeometricContext\geomContext_src_07_02_08\src\geom\');
    vpvar=[0.00001 0.00001 0.00001];
    % chucking out the lines near image boundaries imaging artifacts
    h=imsize(1);
    w=imsize(2);
    inds = find(sum(double(All_lines(:,1:2)>10),2) & sum(double(All_lines(:,1:2)<w-10),2) & ...
        sum(double(All_lines(:,3:4)>10),2) & sum(double(All_lines(:,3:4)<h-10),2));
    All_lines = All_lines(inds,:);
    All_lines=[All_lines sqrt(((All_lines(:,1)-All_lines(:,2)).^2+(All_lines(:,3)-All_lines(:,4)).^2))];
    maxl=max(All_lines(:,7));
      

    %%%%%%Computing intersections of all the lines%%%%%%
    lines = All_lines;
    Xpnts = ComputeIntersectionPoints(lines);
    inds = find(~isnan(Xpnts(:,1)) & ~isnan(Xpnts(:,2)) & ...
        ~isinf(Xpnts(:,1)) & ~isinf(Xpnts(:,2)));
    Xpnts = Xpnts(inds,:);
    VoteArr = ComputeLinePtVote(lines,Xpnts);
    Vote=sum(VoteArr,1);
    [vv ii]=sort(Vote,'descend');
    vp(1:2)=Xpnts(ii(1),1:2);
    Vote1 = VoteArr(:,ii(1));
    active_lines = find((Vote1*maxl./All_lines(:,7))<0.8);
    inactive_lines = find((Vote1*maxl./All_lines(:,7))>=0.8);
    Vote1 = [Vote1(active_lines);Vote1(inactive_lines)];
    lines = All_lines(active_lines,:);

    Xpnts = ComputeIntersectionPoints(lines);
    inds = find(~isnan(Xpnts(:,1)) & ~isnan(Xpnts(:,2)) & ...
        ~isinf(Xpnts(:,1)) & ~isinf(Xpnts(:,2)));
    Xpnts = Xpnts(inds,:);
    VoteArr = ComputeLinePtVote([lines;All_lines(inactive_lines,:)],Xpnts);
    Vote=sum(VoteArr(1:size(lines,1),:),1);
    [vv ii]=sort(Vote,'descend');
    Vote = vv(:);
    Xpnts=Xpnts(ii,:);
    VoteArr = VoteArr(:,ii);
    [Xpnts,Vote,VoteArr] = RemoveRedundantPoints2(Xpnts,Vote,VoteArr,w,h);

    %% Vectorized orthogonality check
    [pts2,pts1]=find(~triu(ones(length(Vote))));
    npts=length(pts1);
    orthochks=[];
    for pt=1:100000:npts
        tempinds = [pt:min(pt+100000-1,npts)];
        temp_orthochks=chckothrogonalityvector(...
            ones(length(tempinds),1)*vp(1:2),...
            Xpnts(pts1(tempinds),:),...
            Xpnts(pts2(tempinds),:),w,h);
        orthochks = [orthochks;temp_orthochks(:)];
    end
    orthos = find(orthochks);
    pts1 = pts1(orthos);
    pts2 = pts2(orthos);
    npts=length(pts1);

    % Total vote computation for these points
    totVote = zeros(npts,1);
    for ln=1:length(Vote1)
    Votes = [Vote1(ln)*ones(npts,1) VoteArr(ln,pts1)' VoteArr(ln,pts2)'];
    Votes = max(Votes,[],2);
    totVote = totVote+Votes;
    end
    totVote = [pts1(:) pts2(:) totVote(:)];

    [vv ii]=sort(totVote(:,3),'descend');
    vp(3:4) = Xpnts(totVote(ii(1),1),:);
    vp(5:6) = Xpnts(totVote(ii(1),2),:);
  
    lines = All_lines;
    VoteArrTemp = ComputeLinePtVote(lines,[vp(1:2);vp(3:4);vp(5:6)]);
    p=[VoteArrTemp.*maxl./repmat(All_lines(:,7),[1 3]) zeros(size(lines,1),1)];%4th vp is outliers
    ind=find(max(p(:,1:3),[],2)< 0.5);
    p(ind,4)=1;
    p=p./repmat(sum(p,2),[1 4]);

    
    [vv linemem] = max(VoteArrTemp,[],2);
    [vv linemem] = max(p,[],2);
    
    %Plot three vps
    if DO_DISPLAY
        figure(1000);plot(vp(1),vp(2),'r*');hold on;
        imagesc(img);hold on;
        plot(vp(1),vp(2),'r*');
        plot(vp(3),vp(4),'g*');
        plot(vp(5),vp(6),'b*');

        % linemem(vv==0) = 4;
        grp1=find(linemem==1);
        grp2=find(linemem==2);
        grp3=find(linemem==3);
        grp4=find(linemem==4);
        plot(lines(grp1, [1 2])', lines(grp1, [3 4])','r');
        plot(lines(grp2, [1 2])', lines(grp2, [3 4])','g');
        plot(lines(grp3, [1 2])', lines(grp3, [3 4])','b');
        plot(lines(grp4, [1 2])', lines(grp4, [3 4])','c');
        axis ij;axis equal;
        saveas(1000,[savedir d(i).name(1:end-3) 'fig']);
        close all

    end

   
        v2(1,:)=vp(1:2);
        v2(2,:)=vp(3:4);
        v2(3,:)=vp(5:6);
        
        v2(:, 1) = (v2(:,1)-imsize(2)/2) / imsize(1);%meansize;
        v2(:, 2) = (v2(:,2) - imsize(1)/2) / imsize(1);%meansize;
        v2(:, 1) = v2(:, 1) + imsize(2)/imsize(1)/2;
        v2(:, 2) = v2(:, 2) + 1/2;
        hpos = APPvp2horizon(vp, vpvar, p, [h w]);
        

    
return