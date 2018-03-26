function imsegs = APPgetSpStats(imsegs)
% imsegs = APPgetSpStats(imsegs)
% Gets basic information about the superpixels
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005

for ii = 1:length(imsegs)
        
	nseg = imsegs(ii).nseg;
	segimage = imsegs(ii).segimage;
	
    imh = size(segimage, 1);
    
	adjmat = eye([nseg nseg], 'uint8');

    % get adjacency
    dx = uint8(segimage ~= segimage(:,[2:end end]));
    dy = segimage ~= segimage([2:end end], :);
            
    ind1 = find(dy);
    ind2 = ind1 + 1;
    s1 = segimage(ind1);
    s2 = segimage(ind2);
    adjmat(s1 + nseg*(s2-1)) = 1;
    adjmat(s2 + nseg*(s1-1)) = 1;
            
    ind3 = find(dx);
    ind4 = ind3 + imh;
    s3 = segimage(ind3);
    s4 = segimage(ind4);
    adjmat(s3 + nseg*(s4-1)) = 1;
    adjmat(s4 + nseg*(s3-1)) = 1;  
    

%   slower code
% 	[height, width] = size(segimage);
% 	
% 	for y = 1:height-1
%         for x = 1:width-1
%             s1 = segimage(y, x);
%             s2 = segimage(y+1, x);
%             s3 = segimage(y, x+1);
%             if s1 > 0
%                 npixels(s1) = npixels(s1) + 1;
%                 if s2 > 0 
%                     adjmat(s1, s2) = 1;            
%                     adjmat(s2, s1) = 1;
%                 end                
%                 if s3 > 0
%                     adjmat(s1, s3) = 1;
%                     adjmat(s3, s1) = 1;
%                 end
%             end
%         end
% 	end
% 	
% 	x = width;
% 	for y = 1:height
%         s1 = segimage(y, x);
%         if s1 > 0
%             npixels(s1) = npixels(s1) + 1;
%         end
% 	end
% 	
% 	y = height;
% 	for x = 1:width-1        
%         s1 = segimage(y, x);
%         if s1 > 0             
%             npixels(s1) = npixels(s1) + 1;
%         end
% 	end

    stats = regionprops(segimage, 'Area');
    imsegs(ii).npixels = vertcat(stats(:).Area);
	imsegs(ii).adjmat = logical(adjmat);
    
end
