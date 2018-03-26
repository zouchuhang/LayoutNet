function [layout]=getLayouts(vp,h,w,nsamp)
% Get candidate box layout. Samples mulitple rays from each vanishing point
%layout is given by corners of the walls.

imcor=[0 h+1;0 0;w+1 0;w+1 h+1];

[vp nothing]=ordervp(vp,h,w);
%inf conditions
infcond = vp(:,1)>50*w | vp(:,2)>50*h;

%inside image conds

inim = vp(:,1)>=1 & vp(:,1)<=w & vp(:,2)>=1 & vp(:,2)<=h;
rays=cell(1,3);


for vno=1:2
    rays{vno}=[];
    
    %general
    if ~inim(vno)
        ll1=[vp(vno,1) imcor(1,1) vp(vno,2) imcor(1,2)];
        ll2=[vp(vno,1) imcor(2,1) vp(vno,2) imcor(2,2)];
        ll3=[vp(vno,1) imcor(3,1) vp(vno,2) imcor(3,2)];
        ll4=[vp(vno,1) imcor(4,1) vp(vno,2) imcor(4,2)];
        lla=[ll1;ll1;ll1;ll2;ll2;ll3];
        llb=[ll2;ll3;ll4;ll3;ll4;ll4];
        
        veca = lla(:,[2,4])-lla(:,[1,3]);
        vecb = llb(:,[2,4])-llb(:,[1,3]);
        norma=(sum(veca.*veca,2)).^.5;
        normb=(sum(vecb.*vecb,2)).^.5;
        theta = acosd(dot(veca,vecb,2)./norma./normb);
        [vv ii]=max(theta);
        
        dtheta=theta(ii)/nsamp;
        rays{vno}(1,:)=lla(ii,:);% two extra rays passing outside the image
        rays{vno}(2,:)=llb(ii,:);
        
        veca=lla(ii,[2,4])-lla(ii,[1,3]);
        vecb=[w 0];
        norma=(sum(veca.*veca,2)).^.5;
        normb=(sum(vecb.*vecb,2)).^.5;
        ang = acosd(dot(veca,vecb,2)./norma./normb);
        
        
        veca=llb(ii,[2,4])-llb(ii,[1,3]);
        vecb=[w 0];
        norma=(sum(veca.*veca,2)).^.5;
        normb=(sum(vecb.*vecb,2)).^.5;
        ang = acosd(dot(veca,vecb,2)./norma./normb);
        
        
    else
        dtheta=180/nsamp;
    end
    
    for ang=dtheta:dtheta:180
        
        pa=[];
        pb=[];
        m=tand(ang);
        
        if isnan(m)
            pa=[vp(vno,1) h];
            pb=[vp(vno,1) 1];
        elseif m==0
            pa=[1 vp(vno,2)];
            pb=[w vp(vno,2)];
        else
            c=vp(vno,2)-m*vp(vno,1);
            x=(h-c)/m;
            p(1,:)=[x h];%intersection with last row of image
            x=(1-c)/m;
            p(2,:)=[x 1];%intersection with 1st row of image
            y=m*w+c;
            p(3,:)=[w y];
            y=m+c;
            p(4,:)=[1 y];
            
            ind=find(p(:,1) <= w & p(:,1)>=1 & p(:,2)<=h & p(:,2)>=1);
            if numel(ind)==2
                pa=p(ind(1),:);
                pb=p(ind(2),:);
            end
        end
        
        if numel(pa) ==2 &numel(pb)==2
            rays{vno}=[rays{vno};pa(1) pb(1) pa(2) pb(2)]; %lines=[x1 x2 y1 y2];
            
        end
    end
end

p1 = [rays{1}(:, [1 3]) ones(size(rays{1}, 1), 1)];
p2 = [rays{1}(:, [2 4]) ones(size(rays{1}, 1), 1)];

ll1 = cross(p1, p2);
ll1 = ll1 ./ repmat(sqrt(sum(ll1.^2,2)), 1, 3);
% ll=cross([vp(2,:) 1],[vp(3,:) 1]);
ll=cross([vp(2,:) 1],[vp(3,:) 1]);
ll = ll ./ repmat(sqrt(sum(ll.^2,2)), 1, 3);
aa = cross(ll1,repmat(ll,[size(ll1,1),1]));
aa=[aa(:,1)./aa(:,3) aa(:,2)./aa(:,3)];
inds = find(aa(:,1)<vp(3,1));
if numel(inds)>0
    rays_left = rays{1}(inds,:);
else
    rays_left = [vp(1,1) vp(3,1)-10 vp(1,2) vp(3,2)];
end
inds = find(aa(:,1)>=vp(3,1));
if numel(inds)>0
    rays_right = rays{1}(inds,:);
else
    rays_right = [vp(1,1) vp(3,1)+10 vp(1,2) vp(3,2)];
end

p1 = [rays{2}(:, [1 3]) ones(size(rays{2}, 1), 1)];
p2 = [rays{2}(:, [2 4]) ones(size(rays{2}, 1), 1)];

ll1 = cross(p1, p2);
ll1 = ll1 ./ repmat(sqrt(sum(ll1.^2,2)), 1, 3);
ll=cross([vp(1,:) 1],[vp(3,:) 1]);
ll = ll ./ repmat(sqrt(sum(ll.^2,2)), 1, 3);
aa = cross(ll1,repmat(ll,[size(ll1,1),1]));
aa=[aa(:,1)./aa(:,3) aa(:,2)./aa(:,3)];
inds = find(aa(:,2)<vp(3,2));
if numel(inds)>0
    rays_top = rays{2}(inds,:);
else
    rays_top = [vp(2,1) vp(3,1) vp(2,2) vp(3,2)-10];
end
inds = find(aa(:,2)>=vp(3,2));
if numel(inds)>0
    rays_bottom = rays{2}(inds,:);
else
    rays_bottom = [vp(2,1) vp(3,1) vp(2,2) vp(3,2)+10];
end

clear rays

[uu vv ww xx]=ndgrid(1:size(rays_left,1),...
    1:size(rays_right,1),...
    1:size(rays_bottom,1),...
    1:size(rays_top,1));
uu=uu(:);
vv=vv(:);
ww=ww(:);
xx=xx(:);


[xs,ys] = IntersectLines(rays_left(uu(:),:),rays_bottom(ww(:),:));
corners_x(:,1) = xs;
corners_y(:,1) = ys;

[xs,ys] = IntersectLines(rays_left(uu(:),:),rays_top(xx(:),:));
corners_x(:,2) = xs;
corners_y(:,2) = ys;

[xs,ys] = IntersectLines(rays_right(vv(:),:),rays_top(xx(:),:));
corners_x(:,3) = xs;
corners_y(:,3) = ys;
[xs,ys] = IntersectLines(rays_right(vv(:),:),rays_bottom(ww(:),:));
corners_x(:,4) = xs;
corners_y(:,4) = ys;


corners_x = round(corners_x);
corners_y = round(corners_y);

ind=find(corners_x >=1 & corners_x <=w & corners_y >=1 & corners_y <=h );
in_img=zeros(size(corners_x));
in_img(ind)=1;

corners_x_temp = corners_x.*in_img;
corners_y_temp = corners_y.*in_img;



[unq,I,J] = unique([corners_x_temp corners_y_temp],'rows');

corners_x = corners_x(I,:);
corners_y = corners_y(I,:);


layout = [];
for i=1:4
    layout = [layout corners_x(:,i) corners_y(:,i)];
end

return;

