function im = displaySegmentsColor(smap, segim, grayim, dodisp, varargin)
% im = displaySegments(smap, segim, grayim, dodisp, options)
% options: 'color': [r g b]
%          'width': line width

if ~exist('dodisp', 'var')
    dodisp = 1;
end

sz=size(segim);
colors = round(255*rand(max(smap),3));
r = zeros(sz);
g = zeros(sz);
b = zeros(sz);
im = zeros(sz(1),sz(2),3);
for i=1:max(smap),
    inds = find(smap(segim)==i);
    r(inds) = round(255*rand(1));
    g(inds) = round(255*rand(1));
    b(inds) = round(255*rand(1));
end
im(:,:,1) = r;
im(:,:,2) = g;
im(:,:,3) = b;


% dispColor = [1 0 0];
% scale = round(max(size(grayim))/200);
% for k = 1:2:numel(varargin)
%     if strcmp(varargin{k}, 'color')
%         dispColor = varargin{k+1};
%     elseif strcmp(varargin{k}, 'width')
%         scale = varargin{k+1};
%     end
% end
% 
%  [h,s,v] = rgb2hsv(label2rgb(smap(segim)));
% % 
%  imagesc(hsv2rgb(h,s,grayim))
% 
% 
% [gx, gy] = gradient(double(smap(segim)));
% g = gx.^2 + gy.^2;
% 
% 
% g = conv2(g, ones(scale), 'same');
% edgepix = find(g>0);
% 
% npix = prod(size(segim));
% 
% im = repmat(grayim, [1 1 3]);
% %im = 0.5*im2double(label2rgb(smap(segim))) + 0.5*im;
% % size(im)
% % size(segim)
% % max(edgepix)
% for b = 1:3
%     im((b-1)*npix+edgepix) = dispColor(b);
% end

if dodisp
    imagesc(im/255), axis image
end