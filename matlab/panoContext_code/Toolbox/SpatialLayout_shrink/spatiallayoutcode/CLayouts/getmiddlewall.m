function [Mwallmask Mwallpolygon Mwalltrapezoid parea]=getmiddlewall(Ra,Rb,Rc,Rd,vp,w,h,flagmask,flagpoly,showpoly)
% get the middle wall polygon given the 4 corners of a box

Mwallmask=[];
Mwallpolygon=[];
Mwalltrapezoid=[];
parea=0;

if nargin < 8
    flagmask=0;
    flagpoly=0;
end

if sum(isnan([Ra Rb Rc Rd]))==0 & flagmask
    Mwallmask=poly2mask([Ra(1) Rb(1) Rc(1) Rd(1)]',[Ra(2) Rb(2) Rc(2) Rd(2)]',h,w);
end

XX1=[Ra(1) Rb(1) Rc(1) Rd(1)]';YY1=[Ra(2) Rb(2) Rc(2) Rd(2)]';


if sum(isnan([Ra Rb Rc Rd]))==0 & flagpoly   %if u need intersection polygon
    XX2=[0 w+1 w+1 0]';YY2=[0 0 h+1 h+1]';
    in=inpolygon(XX1,YY1,[XX2;XX2(1)],[YY2;YY2(1)]);
    
    if numel(find(in==1))==length(in)
        X0=XX1;Y0=YY1;
    else
        [xxx yyy]=meshgrid(1:3);
        xxx=xxx(:);yyy=yyy(:);
        pert=[0 1 -1];
        delx=pert(xxx);delx=delx(:);
        dely=pert(yyy);dely=dely(:);
        
        [in,on]=inpolygon(XX2,YY2,[XX1;XX1(1)],[YY1;YY1(1)]);
        tind=find(on);
        for t=1:numel(tind)
            xv=XX2(tind(t))*ones(size(delx,1),1)+delx;
            yv=YY2(tind(t))*ones(size(dely,1),1)+dely;
            [in,on]=inpolygon(xv,yv,[XX1;XX1(1)],[YY1;YY1(1)]);
            ttind=find(on==0);
            XX2(tind(t))=xv(ttind(1));
            YY2(tind(t))=yv(ttind(1));
        end
        
        
        
        
        [in,on]=inpolygon(XX1,YY1,[XX2;XX2(1)],[YY2;YY2(1)]);;
       
        tind=find(on);
        for t=1:numel(tind)
            xv=XX1(tind(t))*ones(size(delx,1),1)+delx;
            yv=YY1(tind(t))*ones(size(dely,1),1)+dely;
            [in,on]=inpolygon(xv,yv,[XX2;XX2(1)],[YY2;YY2(1)]);
            ttind=find(on==0);
            XX1(tind(t))=xv(ttind(1));
            YY1(tind(t))=yv(ttind(1));
        end
        
        
        [X0 Y0 ind]=polyints(XX1,YY1,XX2,YY2);
    end
    
    if numel(X0)>0 & showpoly
        figure;plot([XX1;XX1(1)],[YY1;YY1(1)],'r');
        hold on;plot([XX2;XX2(1)],[YY2;YY2(1)],'g');
        hold on;plot([X0;X0(1)],[Y0;Y0(1)],'k');
    end
    
    Mwallpolygon=[X0 Y0];
    parea = polyarea(X0,Y0);
    
end

if sum(isnan([Ra Rb Rc Rd]))==0
    Mwalltrapezoid=[XX1 YY1];
end

return;