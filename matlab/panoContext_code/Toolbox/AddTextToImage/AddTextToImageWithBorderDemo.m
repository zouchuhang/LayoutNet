% AddTextToImageWithBorderDemo is a script to demonstrate what
% AddTextToImageWithBorder does.
% It takes about 45 seconds to run on my desktop PC - most of that time is
% spent generating the fonts and displaying the images.
%
% Daniel Warren
% Particle Therapy Cancer Research Institute
% University of Oxford

Image = ones([500 500 3]);

Iterations = 200;

Strings = {'text','image'};

FontNames = {'Arial','Times New Roman','Comic Sans MS','Century Gothic',...
         'Script MT Bold','Courier New','OCR A Extended','Verdana',...
         'Gill Sans MT','Calibri','Cambria'};

FontSizes = 10:20:150;

ThickestBorder = 4;

k = 0;
Fonts = cell(1,length(FontNames)*length(FontSizes));
for i = 1:length(FontNames)
    for j = 1:length(FontSizes)
        k = k+1;
        Fonts{k} = BitmapFont(FontNames{i},FontSizes(j),unique([Strings{:}]));
    end
end

StringIndices = ceil(rand(Iterations,1)*length(Strings));
FontIndices = ceil(rand(Iterations,1)*length(Fonts));
Colors = hsv2rgb([rand(Iterations,1) 0.7+0.3*rand(Iterations,1) 0.9+0.1*rand(Iterations,1)]);
BorderColors = 0.5-(Colors-0.5);
BorderWidths = round(rand(Iterations,1)*ThickestBorder);
XPositions = ceil((-0.2+1.2*rand(1,Iterations))*size(Image,1));
YPositions = ceil((-0.1+1.1*rand(1,Iterations))*size(Image,2));

figure('Units','Pixels','Position',[50 50 50+size(Image,1) 50+size(Image,2)]);
handle = image(Image);
set(gca,'Units','Normalized','Position',[0 0 1 1]);
axis off image;
for i = 1:Iterations
    Image = AddTextToImageWithBorder(Image,Strings{StringIndices(i)},[XPositions(i) YPositions(i)],Colors(i,:),Fonts{FontIndices(i)},[],BorderWidths(i),BorderColors(i,:));
    set(handle,'CData',Image);
    drawnow;
end