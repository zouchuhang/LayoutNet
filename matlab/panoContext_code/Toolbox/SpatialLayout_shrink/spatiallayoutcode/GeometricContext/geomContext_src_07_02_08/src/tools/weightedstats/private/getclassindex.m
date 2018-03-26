function idx = getclassindex(cnames,g)
%GETCLASSINDEX Find indices for class names in another list of names
%   IDX = GETCLASSINDEX(CNAMES,G) takes a list CNAMES of class names
%   (such as the grouping variable values in the treefit or classify
%   function) and another list G of group names (as might be supplied
%   in the "prior" argument to those functions), and finds the indices
%   of the CNAMES names in the G list.  CNAMES should be a cell array
%   of strings.  G can be numbers, a string array, or a cell array of
%   strings

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.1 $  $Date: 2002/03/22 22:03:02 $

% Convert to common string form, whether input is char, cell, or numeric
if isnumeric(g)
   g = cellstr(strjust(num2str(g(:)), 'left'));
elseif ~iscell(g)
   g = cellstr(g);
end

nclasses = length(cnames);
idx = zeros(1,nclasses);

% Look up each class in the grouping variable.
for i = 1:nclasses
   j = strmatch(cnames{i}, g, 'exact');
   if ~isempty(j)
      idx(i) = j(1);
   end
end
