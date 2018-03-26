function [ output ] = smooth( img, sigma )
%SMOOTH Summary of this function goes here
%   Detailed explanation goes here
if isa(img,'uint8')
    img = double(img);
elseif isa(img,'double')
    img = double(img.*255);
end

WIDTH = 4;
sigma = max(sigma, 0.01);
len = 2*ceil(sigma*WIDTH)+1;

padsz = (len-1)/2;
img = padarray(img, [padsz padsz], 'replicate', 'both');

h = fspecial('gaussian', [1 len], sigma);
output = cat(3, conv2(h', h, img(:,:,1),'same'), ...
                conv2(h', h, img(:,:,2),'same'), ...
                conv2(h', h, img(:,:,3),'same'));

output = output(padsz+1:end-padsz, padsz+1:end-padsz, :);
% h = fspecial('gaussian', len, sigma);
% output = filter2(h,img);

end

