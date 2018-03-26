function [ogroup,glabel,gname,multigroup] = mgrp2idx(group,rows,sep);
%MGRP2IDX Convert multiple grouping variables to index vector
%   [OGROUP,GLABEL,GNAME,MULTIGROUP] = MGRP2IDX(GROUP,ROWS) takes
%   the inputs GROUP, ROWS, and SEP.  GROUP is a grouping variable (numeric
%   vector, string matrix, or cell array of strings) or a cell array
%   of grouping variables.  ROWS is the number of observations.
%   SEP is a separator for the grouping variable values.
%
%   The output OGROUP is a vector of group indices.  GLABEL is a cell
%   array of group labels, each label consisting of the values of the
%   various grouping variables separated by the characters in SEP.
%   GNAME is a cell array containing one column per grouping variable
%   and one row for each distinct combination of grouping variable
%   values.  MULTIGROUP is 1 if there are multiple grouping variables
%   or 0 if there are not.

%   Tom Lane, 12-17-99
%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.4 $  $Date: 2002/02/04 19:25:44 $

multigroup = (iscell(group) & size(group,1)==1);
if (~multigroup)
   [ogroup,gname] = grp2idx(group);
   glabel = gname;
else
   % Group according to each distinct combination of grouping variables
   ngrps = size(group,2);
   grpmat = zeros(rows,ngrps);
   namemat = cell(1,ngrps);
   
   % Get integer codes and names for each grouping variable
   for j=1:ngrps
      [g,gn] = grp2idx(group{1,j});
      grpmat(:,j) = g;
      namemat{1,j} = gn;
   end
   
   % Find all unique combinations
   [urows,ui,uj] = unique(grpmat,'rows');
   
   % Create a cell array, one col for each grouping variable value
   % and one row for each observation
   ogroup = uj;
   gname = cell(size(urows));
   for j=1:ngrps
      gn = namemat{1,j};
      gname(:,j) = gn(urows(:,j));
   end
   
   % Create another cell array of multi-line texts to use as labels
   glabel = cell(size(gname,1),1);
   if (nargin > 2)
      nl = sprintf(sep);
   else
      nl = sprintf('\n');
   end
   fmt = sprintf('%%s%s',nl);
   lnl = length(fmt)-3;        % one less than the length of nl
   for j=1:length(glabel)
      gn = sprintf(fmt, gname{j,:});
      gn(end-lnl:end) = [];
      glabel{j,1} = gn;
   end
end
