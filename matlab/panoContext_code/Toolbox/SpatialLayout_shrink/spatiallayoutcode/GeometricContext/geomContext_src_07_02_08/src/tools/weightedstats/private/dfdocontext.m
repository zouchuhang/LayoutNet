function dfdocontext(varargin)
%DFDOCONTEXT Perform context menu actions for distribution fitting tool

% Copyright 2001-2004 The MathWorks, Inc.
% $Revision: 1.1.6.7 $  $Date: 2004/03/09 16:17:04 $
import com.mathworks.toolbox.stats.*;


% Special action to create context menus
if isequal(varargin{1},'create')
   makecontextmenu(varargin{2});
   return
end

% Get information about what invoked this function
obj = gcbo;
action = get(obj,'Tag');
h = gco;
if isempty(h), return; end
dffig = gcbf;

% Set up variables that define some menu items
[sizes styles markers] = getmenuitems;
styles{end+1} = 'none';

changed = true;   % did a line property change?

switch action

 % This case is triggered when we display the menu
 case {'fitcontext' 'datacontext'}
   % Store a handle to the object that triggered this menu
   set(obj,'UserData',h);

   hObject = get(h,'UserData');
   c = findall(obj,'Type','uimenu');
   hBounds = findall(c,'flat','Tag','confbounds');
   ftype = dfgetset('ftype');
   
   % Enable or disable as appropriate
   hMarker = findall(c,'flat','Tag','marker');
   if isequal(action,'datacontext')
      hLineStyle = findall(c,'flat','Tag','linestyle');
      hLineWidth = findall(c,'flat','Tag','linewidth');
      hBinRules = findall(c,'flat','Tag','binrules');
      if isequal(ftype,'probplot')
         set(hMarker,'Enable','on');
         set(hLineStyle,'Enable','off');
         set(hLineWidth,'Enable','off');
         set(hBinRules,'Enable','off');
      elseif isequal(ftype,'pdf')
         set(hMarker,'Enable','off');
         set(hLineStyle,'Enable','on');
         set(hLineWidth,'Enable','on');
         set(hBinRules,'Enable','on');
      else
         set(hMarker,'Enable','off');
         set(hLineStyle,'Enable','on');
         set(hLineWidth,'Enable','on');
         set(hBinRules,'Enable','off');
      end
   else
      if ~hObject.iscontinuous && isequal(ftype,'pdf')
         set(hMarker,'Enable','on');
      else
         set(hMarker,'Enable','off');
      end
   end
   try
      hasconfbounds = hObject.distspec.hasconfbounds;
   catch
      hasconfbounds = false;
   end
   if isequal(ftype,'pdf') || isequal(ftype,'probplot') || ...
      (isequal(action,'fitcontext') && ~hasconfbounds) || ...
      (isequal(action,'datacontext') && isequal(ftype,'icdf'))
      set(hBounds,'Enable','off');
   else
      set(hBounds,'Enable','on');
   end

   set(c,'Checked','off');

   % Fix check mark for confidence bounds
   if hObject.showbounds
      set(hBounds,'Checked','on');
   end

   % Fix check marks on line width and line style cascading menus
   w = get(h,'LineWidth');
   u = findall(c,'flat','Tag',num2str(w));
   if ~isempty(u)
      set(u,'Checked','on');
   end
   w = get(h,'LineStyle');
   u = findall(c,'flat','Tag',w);
   if ~isempty(u)
      set(u,'Checked','on');
   end
   w = get(h,'Marker');
   u = findall(c,'flat','Tag',w);
   if ~isempty(u)
      set(u,'Checked','on');
   end
   return
   
 % Remaining cases are triggered by selecting menu items
 case 'confbounds'
   hObject = get(h,'UserData');
   hObject.showbounds = ~hObject.showbounds;
   nm = get(hObject,'name');
   htag = get(h,'Tag');
   if isequal(htag,'dfdata')
      DataSetsManager.getDataSetsManager.dataSetChanged(java(hObject),nm,nm);
   else
      FitsManager.getFitsManager.fitChanged(java(hObject),nm,nm);
   end

 case 'color'
   oldcolor = get(h,'Color');
   newcolor = uisetcolor(oldcolor);
   if ~isequal(oldcolor,newcolor)
      set(h,'Color',newcolor);
   end

 case styles
   set(h,'LineStyle',action);

 case markers
   if isequal(action,'point')
      msize = 12;
   else
      msize = 6;
   end
   set(h,'Marker',action,'MarkerSize',msize);

 % Either delete a fit, or a hide a fit or data set
 case {'hidecurve' 'deletefit'}
   htag = get(h,'Tag');
   if isequal(htag,'distfit') || isequal(htag,'dfdata')
      hndl = get(h,'UserData');
      if isequal(action,'hidecurve')
         hndl.plot = 0;
         nm = get(hndl,'name');
      else
         % The delete action appears on the fit menu only, not the data set menu
         FitsManager.getFitsManager.deleteFits(java(hndl));
      end
   end
   changed = false;

 % Edit a fit
 case 'editfit'
   htag = get(h,'Tag');
   if isequal(htag,'distfit')  % should always be true
      hndl = get(h,'UserData');
      FitsManager.getFitsManager.editFit(hndl.name);
   end
   changed = false;

 % Bring up the "Set Bin Width Rules" dialog for this data set
 case 'binrules'
   htag = get(h,'Tag');
   if isequal(htag,'dfdata')
      hndl = get(h,'UserData');
      nm = get(hndl,'name');
      bw = com.mathworks.toolbox.stats.BinWidth.getBinWidth; % get dialog
      bw.displayBinWidth(nm); % display dialog for this data set
   end
   changed = false;
 
 % If the menu item is a number, it is a line width
 otherwise
   j = str2num(action);
   if ~isempty(j)
      set(h,'LineWidth',j);
   end

end

if changed
   % Save plot info in the fit or data set object
   hObject = get(h,'UserData');
   savelineproperties(hObject);

   % Update legend
   dfupdatelegend(dffig);
end


% ---------------------- helper to make context menu
function makecontextmenu(dffig)
%MAKECONTEXTMENU Creates context menu for curve fitting figure

% Create context menus for fits, data curve, probability plot data curves
cFit = uicontextmenu('Parent',dffig,'Tag','fitcontext','Callback',@dfdocontext);
uimenu(cFit,'Label','Color...','Tag','color','Callback',@dfdocontext);

% Add menu items for line and marker control
uwidth = uimenu(cFit,'Label','Line &Width','Tag','linewidth');
ustyle = uimenu(cFit,'Label','Line &Style','Tag','linestyle');
umark = uimenu(cFit,'Label','Marker','Tag','marker','Position',2);

% Add menu items to control confidence bounds
uimenu(cFit,'Label','Confidence &Bounds','Callback',@dfdocontext,...
            'Tag','confbounds');

% Get menu item labels and tags
[sizes styles markers slabels mlabels] = getmenuitems;

for j=1:length(markers)
   uimenu(umark,'Label',mlabels{j},'Callback',@dfdocontext,'Tag',markers{j});
end

% Sub-menus for line widths
for i = 1:length(sizes)
   val = num2str(sizes(i));
   uimenu(uwidth,'Label',val,'Callback',@dfdocontext,'Tag',val);
end

% Sub-menus for line styles
for j=1:length(styles)
   uimenu(ustyle,'Label',slabels{j},'Callback',@dfdocontext,'Tag',styles{j});
end

% Copy the fit menu to create a data menu
cData = copyobj(cFit,dffig);
set(cData,'Tag','datacontext')

% Add items for fit menus only
uimenu(cFit,'Label','&Hide Fit','Tag','hidecurve','Callback',@dfdocontext,...
       'Separator','on');
uimenu(cFit,'Label','&Delete Fit','Tag','deletefit','Callback',@dfdocontext);
uimenu(cFit,'Label','&Edit Fit','Tag','editfit','Callback',@dfdocontext);

% Add items for data menus only
uimenu(cData,'Label','&Hide Data','Tag','hidecurve',...
       'Callback',@dfdocontext,'Separator','on');
uimenu(cData,'Label','Set Bin &Rules','Tag','binrules',...
       'Callback',@dfdocontext,'Separator','on');

% -------------- helper to get menu item labels
function [sizes,styles,markers,slabels,mlabels] = getmenuitems
%GETMENUITEMS Get items for curve fitting context menus
sizes = [0.5 1 2 3 4 5 6 7 8 9 10];
styles = {'-' '--' ':' '-.'};
markers = {'+' 'o' '*' '.' 'x' 'square' 'diamond' ...
        'v' '^' '<' '>' 'pentagram' 'hexagram'};
slabels = {'solid' 'dash' 'dot' 'dash-dot'};
mlabels = {'plus' 'circle' 'star' 'point' 'x-mark' 'square' 'diamond' ...
           'triangle (down)' 'triangle (up)' 'triangle (left)' ...
           'triangle (right)' 'pentagram' 'hexagram'};
