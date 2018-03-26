function [sepScene] = separatePano( panoImg, fov, x, y, imgSize, saveDir)
% cut a panorama image into several separate views

if length(x) ~= length(y)
    fprintf('x and y must be the same size.\n');
    return;
end
if length(fov)==1
    fov = fov * ones(length(x),1);
end
% if fov>pi
%     fprintf('I guess the fov is in deg, convert to rad :)\n');
%     fov = fov * pi / 180;
% end
% if ~isdouble('panoImg')
%     fprintf('Image is not double, convert to double and scale to 1 :)');
%     panoImg = double(panoImg)./255;
% end

numScene = length(x);
% imgSize = 2*f*tan(fov/2);
% sepScene = zeros(imgSize, imgSize, 3, numScene);
sepScene(numScene) = struct('img',[],'vx',[],'vy',[],'fov',[],'sz',[]);
parfor i = 1:numScene
    warped_image = imgLookAt(panoImg, x(i), y(i), imgSize, fov(i) );
    sepScene(i).img = warped_image;
    sepScene(i).vx = x(i);
    sepScene(i).vy = y(i);
    sepScene(i).fov = fov(i);
    sepScene(i).sz = imgSize;
end

if exist('saveDir', 'var')
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
        for i = 1:numScene
            imwrite(sepScene(i).img, sprintf('%s\%02d.pgm', saveDir, i), 'pgm');
        end
    end
end