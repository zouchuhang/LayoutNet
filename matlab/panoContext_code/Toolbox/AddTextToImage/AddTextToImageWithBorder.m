function Image = AddTextToImageWithBorder(Image,String,Position,Color,Font,FontSize,BorderWidth,BorderColor)

% Image = AddTextToImage(Image,String,Position,Color,Font,FontSize,BorderWidth,BorderColor)
%
%   Works as AddTextToImage, except also allows for a border of thickness
%   BorderWidth pixels and color given by BorderColor. BorderColor should
%   be a 1- or 3-element vector giving Y or RGB in the range 0 to 1
%   regardless of the image class.
%
%	If no border parameters are provided, default is 1 pixel black.
%
% Daniel Warren
% Particle Therapy Cancer Research Institute
% University of Oxford

if ~exist('BorderWidth','var')
	BorderWidth = 1;
end
if ~exist('BorderColor','var')
	BorderColor = 0;
end

if BorderWidth == 0
	Image = AddTextToImage(Image,String,Position,Color,Font,FontSize);
	return;
end

% uint8 images go from 0 to 255, whereas double ones go from 0 to 1
if isa(Image, 'uint8')
    ScaleFactor = 255;
else
    ScaleFactor = 1;
end

% monochrome images need monochrome text, colour images need colour text
if ndims(Image) == 2 %#ok<ISMAT>
    BorderColor = mean(BorderColor(:));
end
if ndims(Image) == 3 && numel(BorderColor) == 1
    BorderColor = [BorderColor BorderColor BorderColor];
end

Mask = AddTextToImage(false(size(Image(:,:,1))),String,Position,1,Font,FontSize);
ConvElement = double(CircleMask(1+2*BorderWidth,1+2*BorderWidth,BorderWidth+0.5,BorderWidth+0.5,0.5+BorderWidth));
Outline = xor(conv2(double(Mask),ConvElement,'same'),Mask);
Image = AddTextToImage(Image,String,Position,Color,Font,FontSize);

for i = 1:size(Image,3) 
tmp = Image(:,:,i); 
tmp(Outline) = ScaleFactor*BorderColor(i); 
Image(:,:,i) = tmp; 
end

end

function out = CircleMask(w,h,cx,cy,r)
% function out = CircleMask(w,h,cx,cy,r)
%     Outputs a 1D bitmap of a solid circle of radius r pixels
%     centered at pixel (cx,cy) within an matrix of size w x h
%     pixels.

[x y] = meshgrid(1:w,1:h);
out = ((x-cx-0.5).^2+(y-cy-0.5).^2) <= r^2;

end

