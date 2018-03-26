 %Copyright (c) October,15 2008 by Varsha Hedau, UIUC.  All rights reserved.
function   [theta vpinfchk]=LineVpdist(linel,vp)

%%%%%%%%%%%Compute angle between vp and a given line%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%theta in degrees%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%line segment s 
x1=linel(1); x2=linel(2);y1=linel(3);y2=linel(4);
midpntl=[(x1+x2)/2  (y1+y2)/2];
lengthl=sqrt((x1-x2)^2+(y1-y2)^2);

slope1=(y2-y1)/(x2-x1);

%line segment s_dash
slope2=(vp(:,2)-ones(size(vp,1),1)*midpntl(2))./(vp(:,1)-ones(size(vp,1),1)*midpntl(1));
% slope2=(vp(1)-midpntl(1))/(vp(2)-midpntl(2));

%angle between s and s_dash
theta=atand(abs((ones(size(slope2,1),1)*slope1-slope2(:))./(1+ones(size(slope2,1),1)*slope1.*slope2(:))));

%midpoint and slope of s and s_dash are same 
%check if vp lies on rotated s_dash 
d=sqrt((vp(:,1)-ones(size(vp,1),1)*midpntl(1)).^2+(vp(:,2)-ones(size(vp,1),1)*midpntl(2)).^2);
vpinfchk=zeros(size(vp,1),1);
vpinfchk(d > lengthl/2)=1;
return










% function   [theta vpinfchk]=LineVpdist(linel,vp)
% 
% %%%%%%%%%%%Compute angle between vp and a given line%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%theta in degrees%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %line segment s 
% x1=linel(1); x2=linel(2);y1=linel(3);y2=linel(4);
% midpntl=[(x1+x2)/2  (y1+y2)/2];
% lengthl=sqrt((x1-x2)^2+(y1-y2)^2);
% 
% slope1=(y2-y1)/(x2-x1);
% 
% %line segment s_dash
% slope2=(vp(2)-midpntl(2))/(vp(1)-midpntl(1));
% % slope2=(vp(1)-midpntl(1))/(vp(2)-midpntl(2));
% 
% %angle between s and s_dash
% theta=atand(abs((slope1-slope2)/(1+slope1*slope2)));
% 
% %midpoint and slope of s and s_dash are same 
% %check if vp lies on rotated s_dash 
% 
% d=sqrt((vp(1)-midpntl(1))^2+(vp(2)-midpntl(2))^2);
% if d < lengthl/2
% vpinfchk=0;%chk failed if zero dont consider the vote of this segment for the vp 
% else
%     vpinfchk=1;
% end
% return
