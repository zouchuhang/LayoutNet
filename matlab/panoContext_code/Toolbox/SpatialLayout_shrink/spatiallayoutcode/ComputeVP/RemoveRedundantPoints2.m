%Copyright (c) October,15 2008 by Varsha Hedau, UIUC.  All rights reserved.
function [Xpnts2,Vote2,VoteArr2] = RemoveRedundantPoints(Xpnts,Vote,VoteArr,w,h)

%%%%% remove redundant points %%%%%%
currid=1;
done=0;
Xpnts2=Xpnts;
Vote2=Vote;
VoteArr2=VoteArr;

while size(Xpnts,1)>0
    Xpnts2(currid,:) = Xpnts(1,:);
    Vote2(currid) = Vote(1);
    VoteArr2(:,currid) = VoteArr(:,1);
    currid=currid+1;
    dists = (Xpnts(1,1)-Xpnts(:,1)).^2 + ...
        (Xpnts(1,2)-Xpnts(:,2)).^2;
    if sqrt((Xpnts(1,1)-w)^2+(Xpnts(1,2)-h)^2)/sqrt((w/2)^2+(h/2)^2) < 1
        thres=10;
    else
        thres=20*(sqrt((Xpnts(1,1)-w)^2+(Xpnts(1,2)-h)^2)/sqrt((w/2)^2+(h/2)^2));
    end
    inds=find(dists>thres);
    Xpnts = Xpnts(inds,:);
    Vote = Vote(inds,:);
    VoteArr = VoteArr(:,inds);
end
Xpnts2 = Xpnts2(1:currid-1,:);
Vote2 = Vote2(1:currid-1);
VoteArr2 = VoteArr2(:,1:currid-1);


return;
