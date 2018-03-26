function dfadjustmenu(dffig)
%DFADJUSTMENU Adjust contents of curve fit plot menus

%   $Revision: 1.1.6.9 $  $Date: 2004/01/24 09:35:13 $
%   Copyright 2003-2004 The MathWorks, Inc.

% Remove some menus entirely
h = findall(dffig, 'Type','uimenu', 'Parent',dffig);
h0 = findall(h,'flat', 'Label','&Edit');
if (~isempty(h0))
   j = find(h==h0);
   delete(h0);
   h(j) = [];
end
h0 = findall(h,'flat', 'Label','&Insert');
if (~isempty(h0))
   j = find(h==h0);
   delete(h0);
   h(j) = [];
end

% Add or remove some items from other menus
% Fix FILE menu
h0 = findall(h,'flat', 'Label','&File');
h1 = findall(h0, 'Type','uimenu', 'Parent',h0);
m4 = [];
m2 = [];
for j=length(h1):-1:1
   mlabel = get(h1(j),'Label');
   if ~isempty(findstr(mlabel,'Close'))
      m7 = h1(j);
      set(m7,'Label','&Close Distribution Fitting')
   elseif ~isempty(findstr(mlabel,'Print...'))
      m5 = h1(j);
   else
      delete(h1(j));
      h1(j) = [];
   end
end
uimenu(h0, 'Label','&Import Data...', 'Position',1,...
      'Callback','dfittool(''import data'')');
uimenu(h0, 'Label','Clea&r Session','Position',2,...
       'Callback','dfittool(''clear session'')','Separator','on');
uimenu(h0, 'Label','&Load Session...', 'Position',3,...
      'Callback','dfittool(''load session'')');
uimenu(h0, 'Label','&Save Session...', 'Position',4,...
           'Callback','dfittool(''save session'')');
uimenu(h0, 'Label','Generate &M File...', 'Position',5,...
           'Callback','dfittool(''generate code'')');

uimenu(h0, 'Label','&Define Custom Distributions...','Position',6,...
           'Callback',{@dfcustomdist,'define'}','Separator','on');
uimenu(h0, 'Label','I&mport Custom Distributions...', 'Position',7, ...
           'Callback',{@dfcustomdist,'import'},'Tag','importcustom');
uimenu(h0, 'Label','Cl&ear Custom Distributions...', 'Position',8,...
           'Callback',{@dfcustomdist,'clear'},'Tag','clearcustom');


set(m5,'Position',9,'Separator','on');
uimenu(h0, 'Label','Print to &Figure', 'Position',10,...
           'Callback','dfittool(''duplicate'')');
set(m7,'Position',11,'Separator','on');

% Fix VIEW menu
h0 = findall(h,'flat', 'Label','&View');
h1 = findall(h0, 'Type','uimenu', 'Parent',h0);
delete(h1);
uimenu(h0, 'Label','&Legend', 'Position',1,'Separator','off',...
           'Callback','dfittool(''togglelegend'')', 'Checked','on',...
           'Tag','showlegend');
dfgetset('showlegend','on');
uimenu(h0, 'Label','&Grid', 'Position',2,...
           'Callback','dfittool(''togglegrid'')', 'Checked','off', ...
           'Tag','showgrid');
dfgetset('showgrid','off');
h1 = uimenu(h0, 'Label','C&onfidence Level','Position',3,'Separator','on');
uimenu(h1, 'Label','9&0%', 'Position',1, ...
           'Callback','dfittool(''setconflev'',.90)','Tag','conflev');
uimenu(h1, 'Label','9&5%', 'Position',2, 'Checked','on',...
           'Callback','dfittool(''setconflev'',.95)','Tag','conflev');
uimenu(h1, 'Label','9&9%', 'Position',3, ...
           'Callback','dfittool(''setconflev'',.99)','Tag','conflev');
uimenu(h1, 'Label','&Other...', 'Position',4, ...
           'Callback','dfittool(''setconflev'',[])','Tag','conflev');
dfgetset('conflev',0.95);
uimenu(h0, 'Label','&Clear Plot', 'Position',4,...
           'Callback','dfittool(''clear plot'')');

% Fix TOOLS menu
h0 = findall(h,'flat', 'Label','&Tools');
h1 = findall(h0, 'Type','uimenu', 'Parent',h0);
for j=length(h1):-1:1
   mlabel = get(h1(j),'Label');
   if isempty(findstr(mlabel,'Zoom')) && isempty(findstr(mlabel,'Pa&n'))
     delete(h1(j));
     h1(j) = [];
   else
      set(h1(j),'Separator','off');
   end
end
uimenu(h0, 'Label','&Axis Limit Control', 'Position',4, 'Separator','on', ...
           'Callback','dfittool(''toggleaxlimctrl'')', 'Checked','off', ...
           'Tag','showaxlimctrl');
dfgetset('showaxlimctrl','off');
uimenu(h0, 'Label','&Default Axis Limits', 'Position',5, ...
           'Callback','dfittool(''defaultaxes'')');
uimenu(h0, 'Label','Set Default &Bin Rules', 'Position',6, 'Separator','on', 'Callback', ...
           'com.mathworks.toolbox.stats.BinWidth.getBinWidth.displayBinWidth',  ...
           'Tag','setbinrules');
           

% Fix HELP menu
h0 = findall(h,'flat', 'Label','&Help');
h1 = findall(h0, 'Type','uimenu', 'Parent',h0);
delete(h1);
uimenu(h0, 'Label','Statistics &Toolbox Help', 'Position',1,'Callback',...
       'doc stats');
uimenu(h0, 'Label', 'Distribution &Fitting Tool Help', 'Position',2,'Callback',...
        'dfswitchyard(''dfhelpviewer'', ''distribution_fitting'', ''dfittool'')');
uimenu(h0, 'Label','&Demos', 'Position',3,'Separator','on','Callback',...
       'demo toolbox stat'); 
