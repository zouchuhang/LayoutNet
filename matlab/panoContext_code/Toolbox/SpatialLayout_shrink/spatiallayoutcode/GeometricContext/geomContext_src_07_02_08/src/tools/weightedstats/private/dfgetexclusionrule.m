function hExcl = dfgetexclusionrule(ename)
%GETEXCLUSIONRULE Get an exclusion rule by name

% $Revision: 1.1.6.2 $  $Date: 2004/01/24 09:35:37 $
% Copyright 2003-2004 The MathWorks, Inc.

db = getoutlierdb;
hExcl = find(db,'name',ename);

