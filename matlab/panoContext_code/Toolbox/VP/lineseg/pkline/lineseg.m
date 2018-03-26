% LINESEG - Form straight line segements from an edge list.
%
% Usage: [seglist, nedgelist] = lineseg(edgelist, tol, angtol, linkrad)
%
% Arguments:  edgelist - Cell array of edgelists (row col) coords.
%             tol      - Maximum deviation from straight line before a
%                        segment is broken in two (measured in pixels).
%             angtol   - Angle tolerance used when attempting to merge line
%                        segements (radians).
%             linkrad  - Maximum distance between end points of line
%                        segments for segments to be elegible for
%                        linking (pixels).
%  angtol and linkrad are optional.  If these parameters are omitted the
%  merging phase is omitted.
%
% Returns:  seglist - an Nx4 array storing line segments in the form
%                      [x1 y1 x2 y2
%                       x1 y1 x2 y2
%                          . . .    ] etc 
%
%         nedgelist - a new cell array of edge lists where each edge list
%                     corresponds to each segment in seglist above.  The
%                     edgelist is in row,column coords in the form
%                     { [r1 c1   [r1 c1   etc }
%                        r2 c2    ...
%                        ...
%                        rN cN]   ....]
%
% * Note that a non-empty nedgelist is only returned if there is no segment
%   merging. 
%
% This function takes each array of edgepoints in edgelist, finds the
% size and position of the maximum deviation from the line that joins the
% endpoints, if the maximum deviation exceeds the allowable tolerance the
% edge is shortened to the point of maximum deviation and the test is
% repeated.  In this manner each edge is broken down to line segments,
% each of which adhere to the original data with the specified tolerance.
%
% The optional final edge merging phase is provided because the initial
% edge linking phase may have separated edges at `Y' junctions in the
% image.  The merging phase can reconnect broken branches.
% Note however that the merging process can be slow.
%
% See also:  EDGELINK, MAXLINEDEV, MERGESEG, DRAWSEG
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

% December 2000 - Original version
% February 2003 - Added the returning of nedgelist data.

% ** need to preallocate memory to improve speed ***

function [linelist, nedgelist] = lineseg(edgelist, tol, angtol, linkrad)
    
    linelist = [];
    Nline = 0;
    Nedge = length(edgelist);
    
    if nargin == 2
	merge = 0;
    else
	merge = 1;
    end
    
    for e = 1:Nedge
        y = edgelist{e}(:,1);
	x = edgelist{e}(:,2);

	fst = 1;                % Indecies of first and last points in edge
	lst = length(x);        % segment being considered.
	
	while  fst<lst
	    [m,i,d] = maxlinedev(x(fst:lst),y(fst:lst));  % Find size & posn of
                                                          % maximum deviation.
	    
	    while m > tol       % While deviation is > tol  (m/d) ?
		lst = i+fst-1;  % Shorten line to point of max deviation by adjusting lst
		[m,i,d] = maxlinedev(x(fst:lst),y(fst:lst));
	    end
	    
	    Nline = Nline+1;
            % Record line segment. Note that (c,r) corresponds to (x,y)
	    linelist(Nline,:) = [x(fst) y(fst) x(lst) y(lst)]; 
	    % Record edgelist that corresponds to the segment.
	    nedgelist{Nline} = [y(fst:lst) x(fst:lst)];
	               
%	    fst = lst+1;        % Reset fst and lst.
	    fst = lst;        % make new fst match last lst so that segments
                              % share endpoints.
	    lst = length(x);
	end
	
    end
    % fprintf('No of segments = %d\n',length(linelist)); 
    % drawseg(linelist,1), title('Raw segments');
    
    if merge
%	fprintf('\nMerging Segments\n');
	linelist = mergeseg(linelist, angtol, linkrad, tol);
	nedgelist = {};
%	fprintf('No of merged segments = %d\n',length(linelist));
        % drawseg(linelist,2), title('Merged segments');
    end
    

