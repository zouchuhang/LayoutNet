function fitdb=getfitdb(varargin)
% GETFITDB A helper function for DFITTOOL

% $Revision: 1.1.6.3 $  $Date: 2004/01/24 09:36:12 $
% Copyright 2003-2004 The MathWorks, Inc.

thefitdb = dfgetset('thefitdb');

% Create a singleton class instance
if isempty(thefitdb)
   thefitdb = stats.fitdb;
end

dfgetset('thefitdb',thefitdb);
fitdb=thefitdb;
