function [tempimg ]=displayout(Polyg,w,h,img)


r=[255 0 0 0 255 255 0];
g=[0 255 0 255 255 0 0];
b=[0 0 255 255 0 255 0];

rl=[0  0 0 255 255 255 0];
gl=[255 255 0 0 255 0 0];
bl=[0 255 255 0 0 255 0];

fields=zeros(h,w);
%1 floor
polyg=Polyg{1};
if numel(polyg)>0
    tempimg1=poly2mask(polyg(:,1),polyg(:,2),h,w);
else
    tempimg1=zeros(h,w);
end
fields=fields.*(fields~=0)+1*tempimg1.*(fields==0);

%2 middlewall

polyg=Polyg{2};
if numel(polyg)>0
    tempimg1=poly2mask(polyg(:,1),polyg(:,2),h,w);
else
    tempimg1=zeros(h,w);
end
fields=fields.*(fields~=0)+2*tempimg1.*(fields==0);
%3 right wall

polyg=Polyg{3};
if numel(polyg)>0
    tempimg1=poly2mask(polyg(:,1),polyg(:,2),h,w);
else
    tempimg1=zeros(h,w);
end
fields=fields.*(fields~=0)+3*tempimg1.*(fields==0);
%4 left wall

polyg=Polyg{4};
if numel(polyg)>0
    tempimg1=poly2mask(polyg(:,1),polyg(:,2),h,w);
else
    tempimg1=zeros(h,w);
end
fields=fields.*(fields~=0)+4*tempimg1.*(fields==0);
%
%5 ceiling

polyg=Polyg{5};
if numel(polyg)>0
    tempimg1=poly2mask(polyg(:,1),polyg(:,2),h,w);
else
    tempimg1=zeros(h,w);
end
fields=fields.*(fields~=0)+5*tempimg1.*(fields==0);

fields = fields.*(fields~=0) + 6*(fields==0);

mask_r = r(fields);
mask_g = g(fields);
mask_b = b(fields);
mask_color(:,:,1) = mask_r;
mask_color(:,:,2) = mask_g;
mask_color(:,:,3) = mask_b;

tempimg = double(img)*0.5 + mask_color*0.5;
% figure;
% imshow(uint8(tempimg),[]);

return