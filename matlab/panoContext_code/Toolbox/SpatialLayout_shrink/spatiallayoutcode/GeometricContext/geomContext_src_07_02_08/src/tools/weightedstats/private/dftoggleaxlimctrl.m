function dftoggleaxlimctrl(dffig)
%DFTOGGLEAXLIMCTRL Toggle x and y axis limit controls on or off

%   $Revision: 1.1.6.2 $  $Date: 2004/01/24 09:35:53 $
%   Copyright 2003-2004 The MathWorks, Inc.

% Get handle to menu item, may be current object or may not
h = gcbo;
if ~isequal(get(h,'Tag'),'showaxlimctrl')
   h = findall(dffig,'Tag','showaxlimctrl');
end

% Get new state
onoff = on2off(get(h,'Checked'));
dfgetset('showaxlimctrl',onoff);

% Add or remove controls
dfaxlimctrl(dffig,onoff)

% Remove effects of controls on layout
if isequal(onoff,'off')
   dfadjustlayout(dffig);
end

% Change menu state
set(h,'Checked',onoff);

