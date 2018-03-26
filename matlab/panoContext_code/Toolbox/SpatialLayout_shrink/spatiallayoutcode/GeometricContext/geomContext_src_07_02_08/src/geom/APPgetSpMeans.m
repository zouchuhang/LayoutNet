function smeans = APPgetSpMeans(imseg, values)
% smeans = APPgetSpMeans(imseg, values)
% 
% Computes the mean of values for each superpixel in imseg.
%
% imseg.{npixels(nseg), segimage, nseg} - the segmentation information
% values - a matrix of values
% smeans(nseg, 1) - the means of values for each segment

smeans = zeros(imseg.nseg, 1);
count = zeros(imseg.nseg, 1);
segimage = imseg.segimage;
nelements = numel(segimage);
for n = 1:nelements
    s = segimage(n);
    if s~=0
        smeans(s) = smeans(s) + values(n);
        count(s) = count(s) + 1;
    end
end
smeans = smeans ./ (count+ (count==0));


    