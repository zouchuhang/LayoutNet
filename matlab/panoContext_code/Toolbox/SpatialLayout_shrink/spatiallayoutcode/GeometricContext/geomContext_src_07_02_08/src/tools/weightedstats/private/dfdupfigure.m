function dfdupfigure(dffig)
%DFDUPFIGURE Make a duplicate, editable copy of the distribution fitting figure

%   $Revision: 1.1.6.5 $  $Date: 2004/03/26 13:30:58 $
%   Copyright 2003-2004 The MathWorks, Inc.

% Copy the regular axes, not the legend axes
f = figure;
ax = findall(dffig,'Type','axes','Tag','main');
copyobj(ax,f);
newax = findall(f,'Type','axes','Tag','main');

% Adjust layout in new figure, but don't add axis controls
dfadjustlayout(f,'off');

for i=1:length(ax)
   % Remove any context menus and callbacks associated with the old figure
   set(findall(newax(i),'Type','line'),...
       'DeleteFcn','','UIContextMenu',[],'ButtonDownFcn','');

   % Make a new legend based on the original, if any
   [legh,unused,h0,txt] = legend(ax(i));
   if length(h0)>0
      c0 = get(ax(i),'Child');
      c1 = get(newax(i),'Child');
      h1 = h0;
      for j=length(h0):-1:1
         k = find(c0==h0(j));
         if isempty(k)
            h1(j) = [];
            txt(j) = [];
         else
            % Convert to lineseries 
            h1(j) = hgline2lineseries(c1(k(1)));
         end
      end
      
      % It's hard to match the old legend position adequately because the
      % two figure sizes may be very different.  Since the original legend
      % positions were TL or TR, we'll just pick whichever of these is
      % closest to the current position.  If the user moved the legend in
      % the original figure, he or she may need to do the same here.
      oldpos = get(legh,'Position');
      if oldpos(1)<.4
         newpos = 'NW';
      else
         newpos = 'NE';
      end
      legend(newax(i),h1,txt,'Location', newpos);
   end
end
