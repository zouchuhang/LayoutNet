function [ panoout ] = combineViews( Imgs, width, height )
%COMBINEVIEWS Combine separate views to panorama
%   Imgs: same format as separatePano

panoout = zeros(height, width, size(Imgs(1).img,3));
panowei = zeros(height, width, size(Imgs(1).img,3));
imgNum = length(Imgs);
for i = 1:imgNum
    [sphereImg validMap] = ...
        im2Sphere( Imgs(i).img, Imgs(i).fov, width, height, Imgs(i).vx, Imgs(i).vy);
    sphereImg(~validMap) = 0;   
    panoout = panoout + sphereImg;
    panowei = panowei + validMap;
end
panoout(panowei==0) = 0;
panowei(panowei==0) = 1;
panoout = panoout./double(panowei);

end

