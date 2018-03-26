function togglegrid(dffig)
%TOGGLEGGRID Toggle x and y axes grid on or off

%   $Revision: 1.1.6.2 $  $Date: 2004/01/24 09:35:54 $
%   Copyright 2003-2004 The MathWorks, Inc.

% Get new state -- note uimenu state reflects old state, and
% uitoggletool state reflects new state
h = gcbo;
if isequal(get(h,'Type'),'uimenu')
   onoff = on2off(get(h,'Checked'));
else
   onoff = get(h,'State');
end
dfgetset('showgrid',onoff);

% Change grid
ax = findall(dffig,'Type','axes');
for j=1:length(ax)
   if ~isequal(get(ax(j),'Tag'),'legend')
      set(ax(j),'xgrid',onoff,'ygrid',onoff)
   end
end

% Change menu state
h = findall(dffig,'Type','uimenu','Tag','showgrid');
if ~isempty(h), set(h,'Checked',onoff); end

% Change button state
h = findall(dffig,'Type','uitoggletool','Tag','showgrid');
if ~isempty(h), set(h,'State',onoff); end

