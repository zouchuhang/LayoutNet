function [edgeimg lines] = canny_wrapper(img, low, high, minlen, maxerr)

%%
currpath = fileparts(mfilename('fullpath'));

CANNY_BIN = fullfile(currpath, 'Canny/Canny.exe');
TEMPDIR = fullfile(currpath, 'Canny/temp/');
INIMG = fullfile(TEMPDIR, 'in.bmp');
OUTEDGE = fullfile(TEMPDIR, 'edge.txt');
OUTLINE = fullfile(TEMPDIR, 'line.txt');

if ~exist(TEMPDIR,'dir')
    mkdir(TEMPDIR);
end

%%
if nargin < 2, low = 10; end
if nargin < 3, high = 30; end
if nargin < 4, minlen = 20; end
if nargin < 5, maxerr = 1; end

%%
% if length(size(img)) ~= 2
%     error('canny_wrapper: img must be a grayscale image');
% end

%%
imwrite(img, INIMG);

%%
cmd = sprintf('%s %s %s %s %d %d %d %d', ...
    CANNY_BIN, INIMG, OUTEDGE, OUTLINE, low, high, minlen, maxerr);
dos(cmd);

%%
fid = fopen(OUTEDGE, 'rt');
e = textscan(fid, '%d:: x: %d, y: %d');
fclose(fid);

fid = fopen(OUTLINE, 'rt');
l = textscan(fid, 'p1: %d %d, p2: %d %d');
fclose(fid);

%%
e{2} = e{2} + 1;
e{3} = e{3} + 1;
l{1} = l{1} + 1;
l{2} = l{2} + 1;
l{3} = l{3} + 1;
l{4} = l{4} + 1;
%%
edgeimg = zeros(size(img,1),size(img,2));
edgeimg(sub2ind(size(edgeimg), e{3}, e{2})) = 1;

for i=1:length(l{1})
    lines(i).point1 = double([l{1}(i) l{2}(i)]);
    lines(i).point2 = double([l{3}(i) l{4}(i)]);
end
1;
if ~exist('lines','var')
    lines = struct('point1',[],'point2',[]);
end
