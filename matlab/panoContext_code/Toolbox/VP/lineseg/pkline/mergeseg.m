% MERGESEG - Line segment merging function.
%
% Usage: newseglist = mergeseg(seglist, angtol, linkrad)
%
% Arguments:  seglist - an Nx4 array storing line segments in the form
%                       [x1 y1 x2 y2
%                        x1 y1 x2 y2
%                           . . .   ] etc
%
%             angtol  - Angle tolerance used when attempting to merge line
%                       segements (radians).
%             linkrad - Maximum distance between end points of line segments
%                       for segments to be elegible for linking (pixels).
%
% Function scans through the list of line segments seeking to merge segments
% together.  Segments are merged if the orientation difference is less than
% angtol and if the ends of the segments are within linkrad of each other.
% The comparison is performed by first sorting line segments by angle, then
% each line segment is only tested against other line segments that satisfy
% the orientation constraint.
%
% See also:  EDGELINK, LINESEG, MAXLINEDEV, DRAWSEG
%

% Copyright (c) 2000-2005 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% March    2000 - Original version
% November 2005 - Bug fix for case Nseg < 4 (thanks to Ming Luo)
% May      2006 - Bug fix for termination of while loop (thanks to Martin Labrecque)

function newseglist = mergeseg(seglist, angtol, linkrad, linedevtol)
    
    Nseg = size(seglist,1);
    
    if Nseg <= 1    % Nothing to merge
	newseglist = seglist;
	return;
    end
    
    cosAngTol = cos(angtol);
 
    % Build up some data in the following form
    % [x1 y1 x2 y2 dx dy angle merged]
    
    seg = zeros(Nseg,8);
    seg(:,1:4) = seglist;
    seg(:,5:6) = seg(:,3:4)-seg(:,1:2);     % delta x and delta y
    seg(:,7) = atan2(seg(:,6), seg(:,5));   % =segment angle [-pi ... pi]
    neg = seg(:,7) < 0;                     % compress angle range [0...pi]
    seg(:,7) = seg(:,7) + neg*pi;
    [Y,I] = sort(seg(:,7));                 % sort by angle
    seg = seg(I,:);                      
    seg(:,7) = 2*seg(:,7);                  % double the angles
    
    TotalSegs = Nseg;
    
    for s1 = 1:Nseg
%	fprintf('Merging Segment %d/%d\r',s1,Nseg);	    
	if ~seg(s1,8)  % if segment s1 has not already been merged
	    s2 = s1+1; if s2>Nseg, s2=s2-Nseg; end  % `modulo' arithetic
            ang1 = seg(s1,7); ang2 = seg(s2,7);
	    deltaSin = sin(ang1)*cos(ang2) - cos(ang1)*sin(ang2);
	    deltaCos = cos(ang1)*cos(ang2) + sin(ang1)*sin(ang2);
	    deltaAng = 0.5*abs(atan2(deltaSin,deltaCos));

	    while (s2 ~= s1) && (deltaAng < angtol)
		if ~seg(s2,8)   % if segment s2 has not already been merged
		    [newseg, merged] = compareseg(seg(s1,:), seg(s2,:),linkrad, ...
						  linedevtol);
		    if merged
			seg(s1,1:4) = newseg;
			%  update delta x, delta y, and angle ? (I think not)
			seg(s2,8) = merged;  % eliminate s2;
			TotalSegs = TotalSegs - 1;
		    end
		end
		s2 = s2+1; if s2>Nseg, s2=s2-Nseg; end		
		ang2 = seg(s2,7);
		deltaSin = sin(ang1)*cos(ang2) - cos(ang1)*sin(ang2);
		deltaCos = cos(ang1)*cos(ang2) + sin(ang1)*sin(ang2); 
		deltaAng = 0.5*abs(atan2(deltaSin,deltaCos));
	    end
	end
    end

%    fprintf('\n');
    
    % Now clean up list by extracting the unmerged entries
    
    newseglist = zeros(TotalSegs,4);
    newNseg = 0;
    for s = 1:Nseg
	if ~seg(s,8)                 % if this segment has not been merged
	    newNseg = newNseg+1;     % store it in the final list.
	    newseglist(newNseg,:) =  seg(s,1:4);  
	end
    end
    
    
%----------------------------------------------------------------------    
% Internal function that does all the work.  
% Two segments, that already meet the angle constraint, are tested for end
% point proximity.  If they meet this constraint they then tested to see
% what the maximum line deviation would be if they were merged.  If this
% is within tolerance they are merged.
%----------------------------------------------------------------------    
    
function [newseglist, merged] = compareseg(seg1, seg2, linkrad, linedevtol)
    
    s1p1 = seg1(1:2);   % point1 of segment1
    s1p2 = seg1(3:4);   % point2 of segment1
    
    s2p1 = seg2(1:2);   % point1 of segment2
    s2p2 = seg2(3:4);   % point2 of segment2
    
    dir1 = seg1(5:6);
    dir2 = seg2(5:6);
    
    cosAng = dot(dir1,dir2);
    
    if cosAng < 0          % Reverse order of points on seg2 so that both
                           % segments are in the `same' direction
        tmp  = s2p1;
	s2p1 = s2p2;
	s2p2 = tmp;
    end

    merged = 0;                        % Default result: not merged
    newseglist = [seg1; seg2];         % and segments left unchanged
    
    if norm(s1p2-s2p1) < linkrad       % If we satisfy endpoint proximity
	rc = [s1p1; s1p2; s2p1; s2p2];
	[maxdev, index, D] = maxlinedev(rc(:,1),rc(:,2));  % test line deviation
	if maxdev < linedevtol
	    newseglist = [s1p1 s2p2];      % Merge the segments one way
	    merged = 1;
	end
    elseif norm(s2p2-s1p1) < linkrad 
	rc = [s2p1; s2p2; s1p1; s1p2];
	[maxdev, index, D] = maxlinedev(rc(:,1),rc(:,2));
	if maxdev < linedevtol	
	    newseglist = [s2p1 s1p2];      % ... or merge the other way.
	    merged = 1;	
	end
    end	    

    

