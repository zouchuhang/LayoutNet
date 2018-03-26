function [lat,lon] = polyjoin(latcells,loncells)
%POLYJOIN Convert line or polygon parts from cell arrays to vector form
% 
%   [LAT,LON] = POLYJOIN(LATCELLS,LONCELLS) converts polygons from cell
%   array format to column vector format.  In cell array format, each
%   element of the cell array is a vector that defines a separate polygon.
%   A polygon may consist of an outer contour followed by holes separated
%   with NaNs.  In vector format, each vector may contain multiple faces
%   separated by NaNs.  There is no structural distinction between outer
%   contours and holes in vector format.
%
%   See also POLYSPLIT, POLYBOOL, POLYCUT.

% Copyright 1999-2011 The MathWorks, Inc.
% $Revision: 1.4.4.6 $  $Date: 2011/05/17 02:12:10 $

if isempty(latcells) && isempty(loncells)
    lat = reshape([], [0 1]);
    lon = lat;
else
    validateattributes(latcells,{'cell'},{'vector'},mfilename,'LATCELLS',1)
    validateattributes(loncells,{'cell'},{'vector'},mfilename,'LONCELLS',2)
    
    assert(isequal(size(latcells),size(loncells)), ...
        'map:polyjoin:cellvectorSizeMismatch', ...
        '%s and %s must match in size.', ...
        'LATCELLS', 'LONCELLS')
    
    latSizes = cellfun(@size, latcells, 'UniformOutput', false);
    lonSizes = cellfun(@size, loncells, 'UniformOutput', false);
    
    assert(isequal(latSizes,lonSizes), ...
        'map:polyjoin:cellContentSizeMismatch', ...
        'Contents of corresponding cells in %s and %s must match in size.', ...
        'LATCELLS', 'LONCELLS')
    
    M = numel(latcells);
    N = 0;
    for k = 1:M
        N = N + numel(latcells{k});
    end
    
    lat = zeros(N + M - 1, 1);
    lon = zeros(N + M - 1, 1);
    p = 1;
    for k = 1:(M-1)
        q = p + numel(latcells{k});
        lat(p:(q-1)) = latcells{k};
        lon(p:(q-1)) = loncells{k};
        lat(q) = NaN;
        lon(q) = NaN;
        p = q + 1;
    end
    if M > 0
        lat(p:end) = latcells{M};
        lon(p:end) = loncells{M};
    end
end
