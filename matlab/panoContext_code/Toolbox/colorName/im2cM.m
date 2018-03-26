function colorIndex = im2cM(im)

% assume im = double(im);

persistent w2cM;
if isempty(w2cM)
%     load('rectangleDetector/w2cM.mat');
    load('w2cM.mat');
end

index_im = 1+floor(im(:,:,1)/8)+32*floor(im(:,:,2)/8)+32*32*floor(im(:,:,3)/8);

colorIndex = reshape(w2cM(index_im(:)),size(im,1),size(im,2));

end