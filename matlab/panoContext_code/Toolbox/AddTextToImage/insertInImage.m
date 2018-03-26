function imgOut = insertInImage(baseImage,insertionCommand,PVs)
% insertInImage: Embed any text or graphics object in an input image
%
% Use insertInImage to embed, or burn, any text or graphics item into an
% image. You can specify, using a cell array of parameter-value pairs
% (PVs), or using a structure, any valid properties of the specified object
% to insert. ALL ITEMS MUST BE SPECIFIED IN PIXEL COORDINATES, and must not
% extend beyond the edges of the image!
%
% SYNTAX: IMGOUT = INSERTINIMAGE(IMGIN,INSERTIONCOMMAND,PVs)
%
% INPUTS:
%          BASEIMAGE: an image, or a handle to an image (or parent object
%                 containing an image), in which the object is to be
%                 embedded. (The image need not be displayed, unless a
%                 handle is provided.)
%
%          INSERTIONCOMMAND: text, rectangle, line, ellipse, etc. to embed
%                 in the image. Internally, insertInImage calls FEVAL;
%                 anything that works inside an feval command will work
%                 here. For example, you can insert the string 'TESTING' at
%                 [x,y] = [20,30] using feval( @() text('TESTING',20,30]),
%                 so the insertionCommand for this would be:
%                 @() text('TESTING',20,30).
%
%                 TEXT:
%                 @() text(x,y,string)
%
%                 RECTANGLE:
%                 @() rectangle('position',[x y w h])
% 
%                 LINE:
%                 @() line(x,y)
%
%          PVs (OPTIONAL): Cell array or structure of any parameter-value
%                 pairs valid for the TYPE of object you wish to insert.
%                 (Note that this _may_ include a 'position' parameter,
%                 which will overwrite any position set with the insertion
%                 command. For example, when you insert a string, PVs can
%                 be any Parameter-Value pairs valid for TEXT objects. (See
%                 'Text Properties' for details.)
% 
% OUTPUTS:
%          IMGOUT: output RGB image of the same class as imgin, with
%                 embedded text or graphic item(s).
%
% EXAMPLES: 
%
%%% Example 1: Text insertion
%img = imread('rice.png');
%imgOut = insertInImage(img, @()text(100,200,'Test String'),...
%                           {'fontweight','bold','color','r','fontsize',18});
%figure,imshow(imgOut);
%
%%% Example 2: Rectangle insertion
%img = imread('pout.tif');
%f = @() rectangle('position',[55 11 114 120]);
%params = {'linewidth',2,'edgecolor','c'};
%imgOut = insertInImage(img,f,params);
%figure,imshow(imgOut)
%
%%% Example 3: Rectangle, Text, Line
%img = imread('cameraman.tif');
%f = @() rectangle('position',[15 10 155 20]);%[95 145 155 20]
%imgOut = insertInImage(img,f,{'edgecolor','g','linewidth',2});
%f = @() text(20,20,'This is the cameraman')
%imgOut = insertInImage(imgOut,f,{'color','r','fontsize',10,'fontweight','bold'});
%f = @() line([20,40],[35,75]);
%imgOut = insertInImage(imgOut,f,{'color','c','linewidth',2});
%figure,imshow(imgOut)
%
%%%Example 4: Detect, label, and point to yellow circles
%   colors = jet(4);
%   img = imread('coloredChips.png');
%   mask = rgb2gray(img);
%   [centers,radii] = imfindcircles(mask,[20 30],...
%       'Sensitivity',0.89,...
%       'Method','TwoStage',...
%       'ObjectPolarity','Bright');
%   imgOut = insertInImage(img, @()text(20,30,'Auto-detected Yellow Circles'),...
%       {'color',[0.6 0.2 0],'fontsize',24,'fontweight','b'} );%
%   params = {'linewidth',4,'edgecolor','m','linestyle','--','curvature',[1 1]};
%   for ii = 1:size(centers,1)
%       f = @() rectangle('position',...
%           [centers(ii,1)-radii(ii) centers(ii,2)-radii(ii) radii(ii)*2 radii(ii)*2]);
%       imgOut = insertInImage(imgOut,f,params);
%       f = @() line('x',[100,centers(ii,1)],'y',[50,centers(ii,2)]);
%       imgOut = insertInImage(imgOut,f,{'linewidth',3,'color',colors(ii,:)});
%   end
%   figure,imshow(imgOut)
% 
%
% NOTES: 
%          SPECIFYING POSITION: Some items (like TEXT) cannot be created
%          without specifying a position inside the insertionCommand.
%          Others (like RECTANGLE) have default positions, which can be
%          modified using PVs. In the former case, you must include the
%          position with the insertionCommand. In either case, specifying
%          'position' as a PV pair will overwrite the initial position.
%
%          Specifying multiple colors (like EdgeColor and FaceColor) for
%          the same object is not supported.
%
%          This function incorporates and uses functionality I have
%          previously shared on the File Exchage as |createButtonLabel|.
%          
%          Usage of the function requires the Image Processing Toolbox.
% 
%          insertInImage triggers the creation of a temporary figure,
%          which will be momentarily visible during image modification.
%
% Created by Brett Shoelson, Ph.D.
% 10/15/2012
% brett.shoelson@mathworks.com (Comments/suggestions welcome!)
%
% See also: createButtonLabel

%
% Copyright MathWorks, Inc. 2012.

if nargin < 2
    error('INSERTINIMAGE: Requires at least two inputs--a starting image, and an insertion instruction.');
end
if ishandle(baseImage)
    if ~strcmp(get(baseImage,'type'),'image')
        baseImage = findobj(baseImage,'type','image');
        baseImage = baseImage(1);
    end
    baseImage = get(baseImage,'cdata');
end

[m,n,co] = size(baseImage);
if co == 1
    imgOut = repmat(baseImage,[1 1 3]);
elseif co == 3
    imgOut = baseImage;
else
    error('INSERTINIMAGE: Unsupported input image.')
end

% [imgData, bbox, clearMask] = insertObject(insertionCommand,[m,n],PVs);
[imgData, bbox, clearMask] = insertObject(insertionCommand,[m,n]);
imgData = im2double(imgData);
clearMask = logical(clearMask);
for jj = 1:3
    % Zero out anything underlying insertion object (Otherwise, you get semitransparency)
    imgOut(bbox(2):bbox(2)+bbox(4),bbox(1):bbox(1)+bbox(3),jj) = ...
        imgOut(bbox(2):bbox(2)+bbox(4),bbox(1):bbox(1)+bbox(3),jj) .* clearMask;
        % Now paint the image colorplane-by-colorplane
    imgOut(bbox(2):bbox(2)+bbox(4),bbox(1):bbox(1)+bbox(3),jj) = ...
        imgOut(bbox(2):bbox(2)+bbox(4),bbox(1):bbox(1)+bbox(3),jj) + imgData(:,:,jj);
end

switch class(baseImage)
    case 'double'
        imgOut = im2double(imgOut);
    case 'single'
        imgOut = im2single(imgOut);
    case 'logical'
        imgOut = im2bw(imgOut,graythresh(imgOut));
    case 'int16'
        imgOut = im2int16(imgOut);
    case 'uint16'
        imgOut = im2uint16(imgOut);
        %  otherwise: UINT8...no change needed
end

function [imgData, bbox, clearMask] = insertObject(insertionCommand,imSize,PVs)

figure('units','pixels',...
    'tag','tmpcreateButtonLabelfig',...
    'windowstyle','normal',...
    'color',[0 0 0],'Visible','off');
% set(gch, 'Visible', 'off');
% set(findobj('tag','tmpcreateButtonLabelfig'), 'Visible','off');
ax = axes('units','normalized','position',[0 0 1 1],...
    'activepositionproperty','outerposition');
% Create and display matrix of ones same size as input image
workingImg = true(imSize);
imshow(workingImg,'InitialMagnification',100);
hold on
tmp = feval(insertionCommand);
% Apply any input Parameter-Value pairs
% applyPVs(tmp,PVs);
% Temporarily change fontcolor to black for detection, then return it to
% intended color
try
    currColor = get(tmp,'color');
    set(tmp,'color','k');
catch %#ok
    currColor = get(tmp,'edgecolor');
    set(tmp,'edgecolor','k');
end
% Capture region, to determine bounds
F = getframe(ax);
imgData = F.cdata(:,:,1);
imgData = imgData(1:imSize(1), 1:imSize(2));
[y,x] = find(imgData < 1);
try
bbox = [min(x) min(y) max(x)-min(x) max(y)-min(y)] + [-2 -2 4 4];
catch
    fprintf('length x: %d,%d\n', size(x));
    fprintf('length y: %d,%d\n', size(y));
end

imgData = imcrop(imgData,bbox);
% This can be fuzzy, due to the use of getframe. We can create a binary
% representation (and reverse the image as well)...
imgData = imcomplement(imgData);
clearMask = uint8(im2bw(imgData,graythresh(imgData)));
imgData = im2uint8(clearMask*255);
clearMask = 1 - clearMask;
imgData = cat(3,...
    uint8(imgData)*currColor(1),...
    uint8(imgData)*currColor(2),...
    uint8(imgData)*currColor(3));
delete(findobj('tag','tmpcreateButtonLabelfig'));

% NESTED FUNCTIONS
    function applyPVs(obj,pvarray)
        if isempty(pvarray)
            return;
        end
        if isstruct(pvarray)
            set(obj,pvarray);
        else %Cell
            for ii = 1:2:numel(pvarray)
                set(obj,pvarray{ii},pvarray{ii+1});
            end
        end
    end

end

end
