function dfdelgraphexclude
%DFDELGRAPHEXCLUDE Called when the exclusion dataset changes

% Copyright 2003-2004 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2004/01/24 09:35:29 $

% Find the exclusion graph's figure window, and delete it
t = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');
c = get(0,'Child');
f = findobj(c,'flat','Type','figure','Tag','dfexcludegraph');
set(0,'ShowHiddenHandles',t);
delete(f);
