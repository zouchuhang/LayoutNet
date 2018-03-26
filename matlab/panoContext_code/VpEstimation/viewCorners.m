function [ output_args ] = viewCorners( corner, img )
%VIEWCORNERS Summary of this function goes here
%   Detailed explanation goes here
if ~exist('img','var')
    img = zeros(512, 1024);
else
    img = imresize(img, [512 1024]);
    if size(img,3) == 3
        img = 0.7*rgb2gray(img);
    end
end

uv = xyz2uvN(corner);
uvc = uv2coords(uv, 1024, 512);

figure; imshow(img); hold on
color = 'rybg';
for i = 1:size(corner,1)
    plot(uvc(i,1),uvc(i,2),'--ys','LineWidth',0.001,...
                'MarkerFaceColor',color(1),...
                'MarkerSize',4);
%     plot(uvc(i,1),uvc(i,2), ...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',2);
end

end

