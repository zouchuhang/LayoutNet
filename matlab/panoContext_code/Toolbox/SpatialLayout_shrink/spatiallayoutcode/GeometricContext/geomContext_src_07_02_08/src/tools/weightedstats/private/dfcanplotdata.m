function ok = dfcanplotdata(ds,dffig)
%DFCANPLOTDATA Determine if we can plot data in the current plot

%   $Revision: 1.1.6.2 $  $Date: 2004/01/24 09:35:19 $
%   Copyright 2003-2004 The MathWorks, Inc.

% No problem unless we have non-positive data
ok = true;
xlimits = xlim(ds);
xmin = xlimits(1);
if xmin>0
   return
end

% Get handle to control containing the distribution list
if nargin<2
   dffig = dfgetset('dffig');
end
hsel = getappdata(dffig,'selectioncontrols');

% No problem unless we have a probability plot
h = hsel(3);   % handle of display type control
choice = get(hsel(3),'Value');
ftypes = getappdata(hsel(3),'codenames');
ftype = ftypes{choice};
if ~isequal(ftype, 'probplot')
   return
end

% No problem unless distribution has support that excludes some data
ax = get(dffig,'CurrentAxes');
distspec = getappdata(ax,'DistSpec');
lobnd = distspec.support(1);
strict = ~distspec.closedbound(1);
if strict && lobnd>=xmin
   ok = false;
elseif ~strict && lobnd>xmin
   ok = false;
end
