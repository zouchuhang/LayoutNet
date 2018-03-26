function varargout = statgetcolor(ax, linetype, objh)
%STATGETCOLOR Get a color, marker, and linestyle suitable for a new line
%
%   [C,M,L,W] = STATGETCOLOR(AX,LINETYPE,OBJH) gets a color, marker,
%   linestyle, and width for drawing a new line of type LINETYPE (either
%   'data' or 'fit') in the axes AX.  OBJH is the handle for the
%   containing dataset or fit object.

%   Copyright 2003-2004 The MathWorks, Inc.  
%   $Revision: 1.1.6.3 $    $Date: 2004/01/24 09:36:28 $

allcolors = get(ax, 'ColorOrder');
lineproperties = {'Color' 'LineStyle' 'LineWidth' 'Marker'};
allmarkers = {'o' '+' '*' 'x' 's' 'd' '.'};

% For bad call, these will be returned
c = [0 0 0];
m = 'none';
l = '-';
w = 1;

% Get values already stored, if any
if ~isempty(objh) & ~isempty(objh.ColorMarkerLine) ...
                  & nargout<=length(objh.ColorMarkerLine)
   [varargout{1:nargout}] = deal(objh.ColorMarkerLine{1:nargout});
   if nargout>4 && isempty(varargout{5})
      % Supply suitable residual values if they are missing
      varargout(5:nargout) = varargout(1:nargout-4);
      if nargout>=6, varargout{6} = '.'; end
      if nargout>=8, varargout{8} = 1; end
   end
   return
end

switch linetype
 case 'data'
  h = findobj(ax, 'Type','line', 'Tag','dfdata');
  dsdb = dfswitchyard('getdsdb');
  unplottedds = [find(dsdb,'plot',1); find(dsdb,'plot',0)];
  unplottedds(unplottedds==objh) = [];

  % Start data colors from the end, to reduce collisions with fit colors
  ncolors = size(allcolors,1) + 1;
  
  % Find an unused color/marker combination
  for j=1:length(allmarkers)
     m = allmarkers{j};
     h1 = findobj(h,'flat','Marker',m);
     a=1;
     for k=1:size(allcolors,1)
        c = allcolors(ncolors-k,:);
        a = findobj(h1,'flat','Color',c);
        if isempty(a)
           for j=1:length(unplottedds)
              cml = get(unplottedds(j),'ColorMarkerLine');
              if iscell(cml) & ~isempty(cml) & ...
                                isequal(cml{1},c) & isequal(cml{2},m)
                 a = j;
                 break
              end
           end
           if isempty(a)           
              varargout = {c m l w};
              return
           end
        end
     end
  end
  varargout = {c m l w};

 case 'fit'
  w = 2;
  h = findobj(ax, 'Type','line', 'Tag','distfit');
  unplottedfit = [find(getfitdb,'plot',0); find(getfitdb,'plot',1)];

  % Find an unused color/linestyle combination, prefer linestyle = '-'
  allstyles = {'-' '--' '-.'};
  a = 1;
  iter = 0;
  for j=1:length(allstyles)
     l = allstyles{j};
     h1 = findobj(h,'flat','LineStyle',l);
     for k=1:size(allcolors,1)
        iter = iter+1;
        c = allcolors(k,:);
        m = allmarkers{1+mod(iter,length(allmarkers))};
        a = findobj(h1,'flat','Color',c);
        if isempty(a)
           for j=1:length(unplottedfit)
              cml = get(unplottedfit(j),'ColorMarkerLine');
              if iscell(cml) & ~isempty(cml) & ...
                                isequal(cml{1},c) & isequal(cml{3},l)
                 a = j;
                 break
              end
           end
           if isempty(a)           
              varargout = {c m l w};
              return
           end
        end
     end
  end
  
  varargout = {c m l w};
end
