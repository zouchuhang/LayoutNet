function vote = getvote(vp, lines, THRES_THETA)
if length(lines)==0
    vote = 0;
    return;
end

if ~isreal(vp)
    INFDIST = 1e7; % distance to a far away point...
    vp = INFDIST * vp / i;
end

global maxlinelength;
% maxlinelength = max([lines.length]);

% % iterative
% for k = 1:length(lines)
% 	ang = anglebetween(lines(k), vanishpoint);
% 	if ang < VANOPTS.THRES_THETA
% 		vote = vote ...
% 			+ 0.3 * (1 - ang/VANOPTS.THRES_THETA*2) ...
% 			+ 0.7 * (lines(k).length/maxlinelength);
% 	end
% end

% vectorized
point1 = reshape([lines.point1],2,[])';
point2 = reshape([lines.point2],2,[])';
midpoint = (point1 + point2)/2;
% v1 = bsxfun(@minus, vanishpoint, midpoint);
v1 = repmat(vp,size(midpoint,1),1) - midpoint;
% v2 = bsxfun(@minus, point2, midpoint);
v2 = point2 - midpoint;

ang = 180/pi * real(acos( sum(v1.*v2,2) ./ sqrt(sum(v1.^2,2)) ./ sqrt(sum(v2.^2,2)) ));
ang = min([ang 180-ang], [], 2);

linelength = [lines.length]';
vote = sum( (ang<THRES_THETA) .* ...
	(0.3*(1-ang/THRES_THETA) + 0.7*(linelength/maxlinelength)) );
