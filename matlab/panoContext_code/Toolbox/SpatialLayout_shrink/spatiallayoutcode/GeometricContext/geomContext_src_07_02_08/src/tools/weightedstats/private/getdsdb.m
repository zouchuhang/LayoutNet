function thedsdb=getdsdb(varargin)

%   $Revision: 1.1.6.2 $  $Date: 2004/01/24 09:36:11 $
%   Copyright 2003-2004 The MathWorks, Inc.

thedsdb = dfgetset('thedsdb');

% Create a singleton class instance
if isempty(thedsdb)
   thedsdb = stats.dsdb;
   dfgetset('thedsdb',thedsdb);
end


