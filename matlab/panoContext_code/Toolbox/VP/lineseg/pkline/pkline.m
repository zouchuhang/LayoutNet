function lines = pkline(im)
% http://www.csse.uwa.edu.au/~pk/research/matlabfns/

% Find edges using the Canny operator with hysteresis thresholds of 0.1
% and 0.2 with smoothing parameter sigma set to 1.
% edgeim = edge(im,'canny', [0.1 0.2], 1);
edgeim = edge(im,'canny',0.1);
% edgeim = edison_canny(im);

% figure(1), imshow(edgeim); truesize(1)

% Link edge pixels together into lists of sequential edge points, one
% list for each edge contour.  Discard contours less than 10 pixels long.
[edgelist, labelededgeim] = edgelink(edgeim, 30, 8);

% Display the labeled edge image with separate colours for each
% distinct edge (choose your favorite colourmap!)
% figure(2), imagesc(labelededgeim); colormap(vga), 
% axis image, axis off, truesize(2)

% Fit line segments to the edgelists with the following parameters:
tol = 1;         % Line segments are fitted with maximum deviation from
		 % original edge of 2 pixels.
angtol = 0.05;   % Segments differing in angle by less than 0.05 radians
linkrad = 2;     % and end points within 2 pixels will be merged.
[seglist, nedgelist] = lineseg(edgelist, tol, angtol, linkrad);

% Draw the fitted line segments stored in seglist in figure window 3 with
% a linewidth of 2
% drawseg(seglist, 3, 2); axis off

%    figure(1), print -Pjpg -r0 edgeim.jpg
%    figure(2), print -Pjpg -r0 labeledgeim.jpg
%    figure(3), print -Pjpg -r0 segmentim.jpg    

lines = struct('point1', mat2cell(seglist(:,1:2),ones(1,size(seglist,1)),2), ...
	'point2', mat2cell(seglist(:,3:4),ones(1,size(seglist,1)),2));

%% discard short lines
for i = 1:length(lines)
    lines(i).length = norm(lines(i).point1 - lines(i).point2);
end
linelen = [lines.length];
lines = lines(linelen>30);

