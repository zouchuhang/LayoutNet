function [aaa,h] = disp_omap(omap, img, OPACITY)

if ~exist('OPACITY','var')
OPACITY = 0.5;
end

ooo = imresize(double(omap), [size(img,1) size(img,2)], 'nearest');
aaa = uint8(OPACITY*double(img)) + uint8((1-OPACITY)*ooo*256);

% if nargout==0
    h = figure; imshow(aaa)
% elseif nargout==1
% end
