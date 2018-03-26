function Image = AddTextToImage(Image,String,Position,Color,Font,FontSize)

% Image = AddTextToImage(Image,String,Position,Color,Font,FontSize)
%
%   Overlays a rasterized version of the text in String on top of the given
%   Image. The top-left coordinate in pixels is set by Position. Text
%   colour is specified in the variable Color. Font may either be a
%   structure output by BitmapFont or a string specifying a font name.  If
%   the latter, BitmapFont will be called for this font with its size in
%   pixels as specified by FontSize.
%
%   Images may be 1- or 3-channel. Images of class double should have range
%   [0 1] and images of class double should have range [0 255].
%
%   Color specifications should be in the range [0 1] for all RGB and
%   grayscale images regardless of their class.
%
% Daniel Warren
% Particle Therapy Cancer Research Institute
% University of Oxford

if ~exist('Image','var') || isempty(Image)
    % Sample image
    Image = linspace(0,1,500)'*(linspace(0,1,500));
    Image = cat(3,Image,rot90(Image),rot90(Image,2));
end
if ~exist('String','var')
    String = 'No string specified.';
end
if ~exist('Position','var')
    Position = [1 1];
end
if ~exist('Color','var')
    Color = [1 1 0];
end
if ~exist('Font','var')
    Font = 'Arial';
end
if ~exist('FontSize','var')
    FontSize = 32;
end

% uint8 images go from 0 to 255, whereas double ones go from 0 to 1
if isa(Image, 'uint8')
    ScaleFactor = 255;
else
    ScaleFactor = 1;
end

% monochrome images need monochrome text, colour images need colour text
if ndims(Image) == 2 %#ok<ISMAT>
    Color = mean(Color(:));
end
if ndims(Image) == 3 && numel(Color) == 1
    Color = [Color Color Color];
end

% remove overflowing text and/or pad mask to image size

TextMask = RasterizeText(String,Font,FontSize);

% only try adding text if some of it will actually overlay the image

if Position(1) < size(Image,1) && Position(2) < size(Image,2) ...
        && Position(1) + size(TextMask,1) > 0 && Position(2) + size(TextMask,2) > 0

    if Position(1) + size(TextMask,1) > size(Image,1)
        TextMask = TextMask(1:(size(Image,1)-Position(1)),:);
    end
    if Position(2) + size(TextMask,2) > size(Image,2)
        TextMask = TextMask(:,1:(size(Image,2)-Position(2)));
    end
    if any(size(TextMask) ~= [size(Image,1) size(Image,2)]-Position) % save the bottom-right pixel if it's already in the mask
        TextMask(size(Image,1)-Position(1),size(Image,2)-Position(2)) = false;
    end
    
    if Position(1) > 0
        TextMask = cat(1,false(Position(1),size(TextMask,2)),TextMask);
    else
        TextMask = TextMask(1-Position(1):end,:);
    end
    
    if Position(2) > 0
		TextMask = cat(2,false(size(TextMask,1),Position(2)),TextMask);
    else
        TextMask = TextMask(:,1-Position(2):end);
    end
    

    Color = ScaleFactor*Color;

    for i=1:length(Color)
        tmp = Image(:,:,i); % to use logical indexing;
        tmp(TextMask) = Color(i);
        Image(:,:,i) = tmp;
    end

end