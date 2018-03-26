function Font = BitmapFont(Name,Size,Characters,Padding)

% Font = BitmapFont(Name,Size,Characters,Padding)
%
%   Outputs a structure with rasterised binary representations of the
%   non-whitespace characters specified in string Characters. Font used is
%   specified by Name and its size in pixels by Size. Padding is the
%   padding applied after each character in pixels.
%
%   All arguments optional. Default is Courier New at 32 px with 2% kerning
%   with the characters abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXY
%   Z1234567890''�!"�$%&/()=?^�+���,.-<\|;:_>*@#[]{}
%
%   Variable width fonts will work, but fixed width fonts are likely to
%   have fewer kerning issues. Works via screenshots and pops up a figure
%   window, which is unideal and may fail on headless systems: a workaround
%   would be to save pre-generated font files on a desktop machine.
%
%   Will produce better results if font smoothing (ClearType, etc.) is
%   turned off.
%
% Daniel Warren
% Particle Therapy Cancer Research Institute
% University of Oxford

if ~exist('Name','var')
    Name = 'Courier New';
end
if ~exist('Size','var')
    Size = 32;
end
if ~exist('Characters','var')
    Characters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890''�!"�$%&/()=?^�+���,.-<\|;:_>*@#[]{}';
end
if ~exist('Padding','var')
    Padding = 0.02*Size;
end

Size = ceil(Size);
Padding = ceil(Padding);

Bitmaps = cell(1,length(Characters));

% Use a single figure and axis for maximum speed. White background.

fighandle = figure('Position',[50 50 150+Size 150+Size],'Units','pixels','Color',[1 1 1]);
axes('Position',[0 0 1 1],'Units','Normalized');
axis off;

for i = 1:length(Characters)
    % Place each character in the middle of the figure
    texthandle = text(0.5,1,Characters(i),'Units','Normalized','FontName',Name,'FontUnits','pixels','FontSize',Size,'HorizontalAlignment','Center','VerticalAlignment','Top','Interpreter','None','Color',[0 0 0]);
	drawnow;
    % Take a snapshot
    Bitmap = getframe(gcf);
    delete(texthandle);
    % Average RGB to minimise effect of ClearType etc.
    Bitmap = mean(Bitmap.cdata,3);
    % Crop height as appropriate (in MATLAB images, first dimension is
    % height). Some characters will be larger than Size (eg. y and g) -
    % allow for this.
    Bitmap = Bitmap(1:find(mean(Bitmap,2)~=255,1,'last'),:);
    % Crop width to remove all white space
    Bitmap = Bitmap(:,find(mean(Bitmap,1)~=255,1,'first'):find(mean(Bitmap,1)~=255,1,'last'));
    % Pad with kerning value
	Bitmap(:,(end+1):(end+Padding)) = 255;
    % Invert and store in binary format
    Bitmaps{i} = false(size(Bitmap));
    Bitmaps{i}(Bitmap < 160) = true; % This threshold could be changed
end

close(fighandle);

Font.Name = Name;
Font.Size = Size;
Font.Characters = Characters;
Font.Bitmaps = Bitmaps;

end