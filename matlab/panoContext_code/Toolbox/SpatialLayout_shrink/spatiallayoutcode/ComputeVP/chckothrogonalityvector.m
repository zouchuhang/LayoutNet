 %Copyright (c) October,15 2008 by Varsha Hedau, UIUC.  All rights reserved.
%checking orthogonality conditions

%vp=[v1x v1y v2x v2y v3x v3y ]
function [orthochk]=chckothrogonalityvector(vp1s,vp2s,vp3s,w,h)
% vp1s, vp2s, and vp3s multiple vanishing points.

orthochk=zeros(size(vp1s,1),1);

%inf conditions
inf1 = vp1s(:,1)>50*w | vp1s(:,2)>50*h;
inf2 = vp2s(:,1)>50*w | vp2s(:,2)>50*h;
inf3 = vp3s(:,1)>50*w | vp3s(:,2)>50*h;

inds_fff = find(~inf1 & ~inf2 & ~inf3);

inds = find(inf1 & ~inf2 & ~inf3);
temp = vp1s(inds,:);
vp1s(inds,:) = vp3s(inds,:);
vp3s(inds,:) = temp;
inds = find(~inf1 & inf2 & ~inf3);
temp = vp2s(inds,:);
vp2s(inds,:) = vp3s(inds,:);
vp3s(inds,:) = temp;
inds_ffi = find((inf1+inf2+inf3)==1);
% inds_ffi = find(~inf1 & ~inf2 & inf3);

inds = find(inf1 & ~inf2 & inf3);
temp = vp2s(inds,:);
vp2s(inds,:) = vp1s(inds,:);
vp1s(inds,:) = temp;
inds = find(inf1 & inf2 & ~inf3);
temp = vp3s(inds,:);
vp3s(inds,:) = vp1s(inds,:);
vp1s(inds,:) = temp;
inds_fii = find((inf1+inf2+inf3)==2);
% inds_fii = find(~inf1 & inf2 & inf3);

inds_iii = find(inf1 & inf2 & inf3);

if numel(inds_fff)>0
    v1sf = vp1s(inds_fff,:);
    v2sf = vp2s(inds_fff,:);
    v3sf = vp3s(inds_fff,:);
    %     temp_orthochk = ones(size(inds_fff));
    temp_orthochk = zeros(size(inds_fff));
    Mats_11 = v1sf(:,1)+v2sf(:,1);
    Mats_12 = v1sf(:,2)+v2sf(:,2);
    Mats_13 = v1sf(:,1).*v2sf(:,1)+v1sf(:,2).*v2sf(:,2);
    Mats_21 = v1sf(:,1)+v3sf(:,1);
    Mats_22 = v1sf(:,2)+v3sf(:,2);
    Mats_23 = v1sf(:,1).*v3sf(:,1)+v1sf(:,2).*v3sf(:,2);
    Mats_31 = v3sf(:,1)+v2sf(:,1);
    Mats_32 = v3sf(:,2)+v2sf(:,2);
    Mats_33 = v3sf(:,1).*v2sf(:,1)+v3sf(:,2).*v2sf(:,2);

    A_11 = Mats_11-Mats_21; A_12 = Mats_12-Mats_22;
    A_21 = Mats_11-Mats_31; A_22 = Mats_12-Mats_32;
    b_1 = Mats_13-Mats_23; b_2 = Mats_13-Mats_33;
    detA = A_11.*A_22-A_12.*A_21;
    u0 = (A_22.*b_1-A_12.*b_2)./detA;
    v0 = (A_11.*b_2-A_21.*b_1)./detA;
    %     inds = find(u0 > 0.7*w | u0 < 0.3*w | v0 > 0.7*h | v0 < 0.3*h | isnan(u0) | isnan(v0));
%     temp_orthochk(inds) = 0;
   temp = Mats_11.*u0+Mats_12.*v0-Mats_13-u0.*u0-v0.*v0;
    %     inds = find(temp<0 | isnan(temp));
    %     temp_orthochk(inds) = 0;

    f = (temp).^(0.5);

    %     inds = find(f<=100 | f>=5000 | isnan(temp));
    %     temp_orthochk(inds) = 0;
    inds = find(u0 <= 0.7*w & u0 >= 0.3*w & v0 <= 0.7*h & v0 >= 0.3*h & f> 0 & f<=5000 & temp> 0);
    temp_orthochk(inds) = 1;


    %     temp_orthochk2 = zeros(size(temp_orthochk));
    %     for i=1:length(temp_orthochk2)
    %         temp_orthochk2(i) = chckothrogonality([v1sf(i,:) v2sf(i,:) v3sf(i,:)],w,h);
    %     end

    orthochk(inds_fff) = temp_orthochk;

end

if numel(inds_ffi)>0
    v1sf = vp1s(inds_ffi,:);
    v2sf = vp2s(inds_ffi,:);
    v3si = vp3s(inds_ffi,:);
%     temp_orthochk = ones(size(inds_ffi));
    temp_orthochk = zeros(size(inds_ffi));
    r=((w/2-v1sf(:,1)).*(v2sf(:,1)-v1sf(:,1))+(h/2-v1sf(:,2)).*(v2sf(:,2)-v1sf(:,2)))./((v2sf(:,1)-v1sf(:,1)).^2+(v2sf(:,2)-v1sf(:,2)).^2);
%     
%     inds = find(r<=0 | r>=1 | isnan(r));
%     temp_orthochk(inds) = 0;

    u0= v1sf(:,1) + r.*(v2sf(:,1)-v1sf(:,1));
    v0= v1sf(:,2) + r.*(v2sf(:,2)-v1sf(:,2));
%     inds = find(u0 > 0.7*w | u0 < 0.3*w | v0 > 0.7*h | v0 < 0.3*h | isnan(u0) | isnan(v0));
%     temp_orthochk(inds) = 0;

    temp=u0.*(v1sf(:,1)+v2sf(:,1))+v0.*(v2sf(:,2)+v1sf(:,2))-(v1sf(:,1).*v2sf(:,1)+v2sf(:,2).*v1sf(:,2)+u0.^2+v0.^2);
%     inds = find(temp<0 | isnan(temp));
%     temp_orthochk(inds) = 0;
    f = (temp).^(0.5);
%     inds = find(f<=0 | f>=5000 | isnan(f));
%     temp_orthochk(inds) = 0;

    vec1=[v2sf(:,1)-v1sf(:,1) v2sf(:,2)-v1sf(:,2)]; vec2=[v3si(:,1) v3si(:,2)];
    dot12 = sum(vec1.*vec2,2);
    norm1 = (sum(vec1.*vec1,2).^.5);
    norm2 = (sum(vec2.*vec2,2).^.5);
%     inds = find(abs(dot12./norm1./norm2)>=0.1 | isnan(dot12) | isnan(norm1) | isnan(norm2));
%     temp_orthochk(inds) = 0;


 inds=find(r > 0 & r < 1 & u0 <= 0.7*w & u0 >= 0.3*w & v0 <= 0.7*h & v0 >= 0.3*h & temp > 0 & f> 0 & f<=5000 ...
     & abs(dot12./norm1./norm2)< 0.1);
    temp_orthochk(inds) = 1;

    %     temp_orthochk2 = zeros(size(temp_orthochk));
    %     for i=1:length(temp_orthochk2)
    %         temp_orthochk2(i) = chckothrogonality([v1sf(i,:) v2sf(i,:) v3sf(i,:)],w,h);
    %     end

    orthochk(inds_ffi) = temp_orthochk;
end

if numel(inds_fii)>0
    v1sf = vp1s(inds_fii,:);
    v2si = vp2s(inds_fii,:);
    v3si = vp3s(inds_fii,:);
%     temp_orthochk = ones(size(inds_fii));

    temp_orthochk = zeros(size(inds_fii));

    vec1 = [v2si(:,1) v2si(:,2)];
    vec2 = [v3si(:,1) v3si(:,2)];
    dot12 = sum(vec1.*vec2,2);
    norm1 = (sum(vec1.*vec1,2).^.5);
    norm2 = (sum(vec2.*vec2,2).^.5);
    
%     inds = find(abs(dot12./norm1./norm2)>=0.1 | isnan(dot12) | isnan(norm1) | isnan(norm2));
%     temp_orthochk(inds) = 0;
% 
%     inds = find(v1sf(:,1)<0.3*w | v1sf(:,1)>0.7*w | v1sf(:,2)<0.3*h | v1sf(:,2)>0.7*h | isnan(v1sf(:,1)) | isnan(v1sf(:,2)));
%     temp_orthochk(inds) = 0;

inds = find(abs(dot12./norm1./norm2)<=0.1 & v1sf(:,1)> 0 & v1sf(:,1)<=0.7*w & v1sf(:,2)> 0 & v1sf(:,2)<=0.7*h );
temp_orthochk(inds) = 1;


    %     temp_orthochk2 = zeros(size(temp_orthochk));
    %     for i=1:length(temp_orthochk2)
    %         temp_orthochk2(i) = chckothrogonality([v1sf(i,:) v2sf(i,:) v3sf(i,:)],w,h);
    %     end

    orthochk(inds_fii) = temp_orthochk;
end

if numel(inds_iii)>0
    orthochk(inds_iii) = 0;
end

