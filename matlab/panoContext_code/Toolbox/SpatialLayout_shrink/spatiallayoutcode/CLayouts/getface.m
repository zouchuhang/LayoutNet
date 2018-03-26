function [face_mask face_polygon parea]=getface(R1,R2,vp,w,h,flagmask,flagpoly,showpoly)

face_mask = [];
face_polygon = [];
parea = 0;

vp = vp(3,:);

polygon = [];

% Intersect the ray 'vp-R1' with the image
[lineout] = LineEndPtsInImg([vp(1) R1(1) vp(2) R1(2)],w,h);
ptsout = reshape(lineout,2,2);
ptsout=round(ptsout);
% first select the points which are on the image bndy and are valid
inds = find((ptsout(:,1)>=1 & ptsout(:,2)>=1 & ...
    ptsout(:,1)<=w & ptsout(:,2)<=h & ...
    ~isnan(ptsout(:,1)) & ~isnan(ptsout(:,2)) & ...
    ~isinf(ptsout(:,1)) & ~isinf(ptsout(:,2))));
if numel(inds)>0
    ptsout = ptsout(inds,:);
    % then select the points which are the same side of vp as R1
    norma = (sum((ptsout-repmat(vp,size(ptsout,1),1)).^2,2)).^.5;
    normb = (sum((R1-vp).^2)).^.5;
    inds = find((sign(ptsout(:,1)-vp(1))==sign(R1(1)-vp(1))) & ...
        (sign(ptsout(:,2)-vp(2))==sign(R1(2)-vp(2))) & (norma>=normb));
    if numel(inds)>0
        polygon = [polygon;ptsout(inds,:)];
    end
end

% Intersect the ray 'vp-R2' with the image
[lineout] = LineEndPtsInImg([vp(1) R2(1) vp(2) R2(2)],w,h);
ptsout = reshape(lineout,2,2);
ptsout=round(ptsout);
% first select the points which are on the image bndy and are valid
inds = find((ptsout(:,1)>=1 & ptsout(:,2)>=1 & ...
    ptsout(:,1)<=w & ptsout(:,2)<=h & ...
    ~isnan(ptsout(:,1)) & ~isnan(ptsout(:,2)) & ...
    ~isinf(ptsout(:,1)) & ~isinf(ptsout(:,2))));
if numel(inds)>0
    ptsout = ptsout(inds,:);
    % then select the points which are the same side of vp as R1
    norma = (sum((ptsout-repmat(vp,size(ptsout,1),1)).^2,2)).^.5;
    normb = (sum((R2-vp).^2)).^.5;
    inds = find((sign(ptsout(:,1)-vp(1))==sign(R2(1)-vp(1))) & ...
        (sign(ptsout(:,2)-vp(2))==sign(R2(2)-vp(2))) & (norma>normb));
    if numel(inds)>0
        polygon = [polygon;ptsout(inds,:)];
    end
end

% After getting the points on the image boundary, select among the image
% corners that should belong to the polygon
ray1 = [vp(1) R1(1) vp(2) R1(2)];
ray2 = [vp(1) R2(1) vp(2) R2(2)];
imcor=[1 h;1 1;w 1;w h];

for i=1:size(imcor,1)
    ray = [vp(1) imcor(i,1) vp(2) imcor(i,2)];
    
    veca = ray1([2,4])-ray1([1,3]);
    vecb = ray([2,4])-ray([1,3]);
    norma=(sum(veca.*veca,2)).^.5;
    normb=(sum(vecb.*vecb,2)).^.5;
    theta1 = acosd(dot(veca,vecb,2)./norma./normb);
    
    veca = ray2([2,4])-ray2([1,3]);
    vecb = ray([2,4])-ray([1,3]);
    norma=(sum(veca.*veca,2)).^.5;
    normb=(sum(vecb.*vecb,2)).^.5;
    theta2 = acosd(dot(veca,vecb,2)./norma./normb);
    
    veca = ray1([2,4])-ray1([1,3]);
    vecb = ray2([2,4])-ray2([1,3]);
    norma=(sum(veca.*veca,2)).^.5;
    normb=(sum(vecb.*vecb,2)).^.5;
    theta3 = acosd(dot(veca,vecb,2)./norma./normb);
    cond1 = abs(theta3-theta2-theta1)<0.01;
    
    det1 = det([R1 1;R2 1;vp 1]);
    det2 = det([R1 1;R2 1;imcor(i,:) 1]);
    cond2 = sign(det1)~=sign(det2);
    
    if cond1 & cond2
        polygon = [polygon;imcor(i,:)];
    end
    
end

if numel(polygon)>0

    polygon = [polygon;R1];
    polygon = [polygon;R2];
    
    [unq,I,J] = unique([polygon(:,1) polygon(:,2)],'rows');
    polygon = polygon(I,:);
    
    
    
    
    if size(polygon,1) > 2 & ~chckcollinearity(polygon)
        K = convhull(polygon(:,1),polygon(:,2));
        XX1 = [polygon(K,1)];
        XX1=XX1(:);
        YY1 = [polygon(K,2)];
        YY1=YY1(:);
    else
       return; 
    end
    
   


    if exist('XX1','var') & sum(isnan([XX1;YY1]))==0  & flagpoly   %if u need intersection polygon

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

            if numel(X0)>0 & showpoly
                figure;plot([XX1;XX1(1)],[YY1;YY1(1)],'r');
                hold on;plot([XX2;XX2(1)],[YY2;YY2(1)],'g')
                hold on;plot([X0;X0(1)],[Y0;Y0(1)],'k');
                axis equal
            end
       
        end
        face_polygon(:,1)=X0;
        face_polygon(:,2)=Y0;
        parea = polyarea(face_polygon(:,1),face_polygon(:,2));
    end
end
return;
