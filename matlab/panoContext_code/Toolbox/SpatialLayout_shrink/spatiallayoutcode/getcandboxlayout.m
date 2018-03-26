function [ polyg, Features] = getcandboxlayout( vp,h,w, integData )
% getcandoxlayout Get box layout candidates and their features.
%
% INPUT:
%    vp- vanishing points h,w- height width of image
%    integData - integral images for features.

%OUTPUT:
% polyg - box layout candidates, each candidate as five polygons corresponding to left/middle/right walls, floor and ceiling
% Features -for each box layout candidate

%Features:-
%line-membership features:1:20
%correct label confidence: 21-25
%label entropy :26 -30
%weighted line-memership features:31-40
%surface label confidence :41 -75


numf=75;

samps=10;gtexist=0;
[layout]=getLayouts(vp,h,w,samps);

numL=size(layout,1);

Features=zeros(numL,numf);
polyg=cell(numL,5);

cnt=1;quantsiz=500;
if numL > 0
    for lay=1:size(layout,1)
        Ra=[layout(lay,1) layout(lay,2)];
        Rb=[layout(lay,3) layout(lay,4)];
        Rc=[layout(lay,5) layout(lay,6)];
        Rd=[layout(lay,7) layout(lay,8)];
        Polyg=[];
        tot_area = 0;
        
        cond1=inpolygon(vp(3,1),vp(3,2),[Ra(1); Rb(1); Rc(1); Rd(1); Ra(1)],[Ra(2); Rb(2); Rc(2); Rd(2);Ra(2)]) & ...
            ~(inpolygon(vp(2,1),vp(2,2),[Ra(1); Rb(1); Rc(1); Rd(1);Ra(1)],[Ra(2); Rb(2); Rc(2); Rd(2);Ra(2)])) & ...
            ~(inpolygon(vp(1,2),vp(1,2),[Ra(1); Rb(1); Rc(1); Rd(1);Ra(1)],[Ra(2); Rb(2); Rc(2); Rd(2);Ra(2)]));
        
        lines1=[Ra(1) Rb(1) Ra(2) Rb(2);...
            Rc(1) Rd(1) Rc(2) Rd(2);...
            Rb(1) Rc(1) Rb(2) Rc(2);...
            Ra(1) Rd(1) Ra(2) Rd(2)];
        lines2=[vp(2,1) vp(3,1) vp(2,2) vp(3,2);...
            [vp(2,1) vp(3,1) vp(2,2) vp(3,2)];...
            [vp(1,1) vp(3,1) vp(1,2) vp(3,2)];...
            [vp(1,1) vp(3,1) vp(1,2) vp(3,2)]];
        [inxs inys]=IntersectLines(lines1,lines2);
        
        cond2=inxs(1) < vp(3,1) & inxs(2) > vp(3,1) & inys(3) < vp(3,2) & inys(4) > vp(3,2);
        
        %1 floor
        [tempimg tpolyg parea]=getface(Rd,Ra,vp,w,h,0,1,0);
        Polyg{1}=[tpolyg];
        if numel(find(isnan(tpolyg))) >0
            error('debug me');
        end
        tot_area = tot_area+parea;
        
        %2 middlewall
        [tempimg tpolyg trapez parea]=getmiddlewall(Ra,Rb,Rc,Rd,vp,w,h,0,1,0);
        Polyg{2}=[tpolyg];
        if numel(find(isnan(tpolyg))) >0
            error('debug me');
        end
        tot_area = tot_area+parea;
        
        %3 right wall
        [tempimg tpolyg parea]=getface(Rc,Rd,vp,w,h,0,1,0);
        Polyg{3}=[tpolyg];
        if numel(find(isnan(tpolyg))) >0
            error('debug me');
        end
        tot_area = tot_area+parea;
        
        
        %4 left wall
        [tempimg tpolyg parea]=getface(Ra,Rb,vp,w,h,0,1,0);
        Polyg{4}=[tpolyg];
        if numel(find(isnan(tpolyg))) >0
            error('debug me');
        end
        tot_area = tot_area+parea;
        
        %5 ceiling
        [tempimg tpolyg parea]=getface(Rb,Rc,vp,w,h,0,1,0);
        Polyg{5}=[tpolyg];
        if numel(find(isnan(tpolyg))) >0
            error('debug me');
        end
        tot_area = tot_area+parea;
        cond3=(numel(Polyg{1})+numel(Polyg{2})+numel(Polyg{3})+numel(Polyg{4})+numel(Polyg{5})) > 0;
        
        
        if (tot_area/w/h)<0.8 | (tot_area/w/h) > 1.2 | ~cond1 | ~cond2 | ~cond3
            continue;
        else
            
            for f=1:5
                polyg{cnt,f}=Polyg{f};
            end
            [Polyg]=clipPolyg(Polyg,h,w);
            if numel(Polyg)>0
                %                     tic
                [features] =getLayoutfeats(Polyg,integData,vp,h,w,0,quantsiz);
                %                     toc
            end
            
            Features(cnt,:)=features;
            
            cnt=cnt+1;
        end
        
        
    end
    
    Features(cnt:end,:)=[];
    polyg(cnt:end,:)=[];
    
    if size(Features,1)~=size(polyg,1)
        disp('stop');
        keyboard;
    end
    
    
end
end






