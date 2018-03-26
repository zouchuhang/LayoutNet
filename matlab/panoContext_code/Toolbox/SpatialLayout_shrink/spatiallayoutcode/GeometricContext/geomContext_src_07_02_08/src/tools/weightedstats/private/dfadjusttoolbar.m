function dfadjusttoolbar(dffig)
%DFADJUSTTOOLBAR Adjust contents of distribution fitting plot toolbar

%   $Revision: 1.1.6.5 $  $Date: 2004/01/24 09:35:14 $
%   Copyright 2003-2004 The MathWorks, Inc.

h0 = findall(dffig,'Type','uitoolbar');
h1 = findall(h0,'Parent',h0);
czoom = [];
for j=length(h1):-1:1
   mlabel = get(h1(j),xlate('TooltipString'));
   if ~isempty(findstr(mlabel,'Zoom')) || ~isempty(findstr(mlabel,'Pan'))
      czoom(end+1) = h1(j);
   elseif isempty(findstr(mlabel,'Print'))
      delete(h1(j));
      h1(j) = [];
   else
     c1 = h1(j);
   end
end

% Add more icons especially for distribution fitting
if exist('dficons.mat','file')==2
   icons = load('dficons.mat','icons');
   state = dfgetset('showlegend');
   if isempty(state), state = 'on'; end
   try
      % Try to get the default MATLAB legend icon
      legicon = load([matlabroot '/toolbox/matlab/icons/legend.mat']);
      cdata = legicon.cdata;
   catch
      cdata = icons.icons.legend; % in case of trouble, use older icon
   end
   c2 = uitoggletool(h0, 'CData',cdata,...
                    'State',state,...
                    'TooltipString', 'Legend On/Off',...
                    'Separator','on',...
                    'ClickedCallback','dfittool(''togglelegend'')',...
                    'Tag','showlegend');
   state = dfgetset('showgrid');
   if isempty(state), state = 'off'; end
   c3 = uitoggletool(h0, 'CData',icons.icons.grid,...
                    'State',state,...
                    'TooltipString', ('Grid On/Off'),...
                    'Separator','off',...
                    'ClickedCallback','dfittool(''togglegrid'')',...
                    'Tag','showgrid');
   c = get(h0,'Children');
   cnew = [c1 czoom c2 c3]';
   
   set(h0,'Children',cnew(end:-1:1));
end
