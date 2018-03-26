function toggletoolbar(varargin)
%TOGGLETOOLBAR Toggle distribution fit plot toolbar on or off

%   $Revision: 1.1.6.2 $  $Date: 2004/01/24 09:35:56 $
%   Copyright 2003-2004 The MathWorks, Inc.

if (nargin>0 & ishandle(varargin{1}) & ...
               isequal(get(varargin{1},'Type'),'figure'))
   dffig = varargin{1};
else
   dffig = gcbf;
end

tbstate = get(dffig,'toolbar');
h = findall(dffig,'Type','uitoolbar');
if isequal(tbstate,'none') | isempty(h)
   % Create toolbar for the first time
   set(dffig,'toolbar','figure');
   dfadjusttoolbar(dffig);
elseif nargin>1 & isequal(varargin{2},'on')
   % Hide toolbar
   set(h,'Visible','on');
else
   % Show toolbar
   set(h,'Visible','off');
end
