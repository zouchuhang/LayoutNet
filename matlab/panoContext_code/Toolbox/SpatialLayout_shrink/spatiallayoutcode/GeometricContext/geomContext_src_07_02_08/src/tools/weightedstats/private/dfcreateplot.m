function dffig = dfcreateplot
%DFCREATEPLOT Create plot window for DFITTOOL

%   $Revision: 1.1.6.8 $  $Date: 2004/03/09 16:17:03 $
%   Copyright 2003-2004 The MathWorks, Inc.

% Get some screen and figure position measurements
tempFigure=figure('visible','off','units','pixels',...
                  'Tag','Distribution Fitting Figure');
dfp=get(tempFigure,'position');
dfop=get(tempFigure,'outerposition');
diffp = dfop - dfp;
xmargin = diffp(3);
ymargin = diffp(4);
close(tempFigure)
oldu = get(0,'units');
set(0,'units','pixels');
screenSize=get(0,'screensize');
screenWidth=screenSize(3); 
screenHeight=screenSize(4);
set(0,'units',oldu');

% Get the desired width and height
width=dfp(3)*1.2 + xmargin;
height=dfp(4)*1.2 + ymargin;
if width > screenWidth
  width = screenWidth-10-xmargin;
end
if height > screenHeight;
  height = screenHeight-10-ymargin;
end

% Calculate the position on the screen
leftEdge=min((screenWidth/3)+10+xmargin/2, screenWidth-width-10-2*xmargin);
bottomEdge=(screenHeight-height)/2;

% Make an invisible figure to start
dffig=figure('Visible','off','IntegerHandle','off',...
             'HandleVisibility','callback',...
             'color',get(0,'defaultuicontrolbackgroundcolor'),...
             'name','Distribution Fitting Tool',...
             'numbertitle','off',...
             'units','pixels',...
             'position',[leftEdge bottomEdge width height], ...
             'CloseRequestFcn',@closefig, ...
             'PaperPositionMode','auto',...
             'doublebuffer','on',...
             'Dock','off');

dfgetset('dffig',dffig);

% Set default print options
pt = printtemplate;
pt.PrintUI = 0;
set(dffig,'PrintTemplate',pt)

% Add buttons along the top
dfaddbuttons(dffig);

% We want a subset of the usual toolbar
% Instead of calling adjusttoolbar, there is a handlegraphics bug
% that turned the toolbar off when the buttons were created, so
% we have to toggle it back on.
dftoggletoolbar(dffig,'on');

% We want a subset of the usual menus and some more toolbar icons
dfadjustmenu(dffig);

% Set up axes the way we want them
ax=axes('Parent',dffig, 'box','on','Tag','main',...
        'XLimMode','manual','YLimMode','manual','ZLimMode','manual',...
        'CLimMode','manual','AlimMode','manual');

% Adjust layout of buttons and graph
if ~ispc    % some unix platforms seem to require this
   set(dffig,'Visible','on');
   drawnow;
end
dfadjustlayout(dffig);

% Remember current position
dfgetset('oldposition',get(dffig,'Position'));
dfgetset('oldunits',get(dffig,'Units'));

% Now make the figure visible
if ispc
   set(dffig,'visible','on');
end
set(dffig, 'ResizeFcn','dfittool(''adjustlayout'')');
drawnow;

% Set up some listeners
hgpkg = findpackage('hg');
axesC = hgpkg.findclass('axes');

% Create context menus for data and fit lines
dfdocontext('create', dffig);

% Listen for figure position changes if resize function is questionable
if ~ispc
   list(1) = handle.listener(dffig, findprop(handle(dffig),'position'), ...
             'PropertyPostSet', 'dfittool(''adjustlayout2'')');
   dfgetset('figlistener',list);
end


% ---------------------- helper to verify closing of figure
function closefig(varargin)
%CLOSEFIG Verify intention to close distribution fitting figure

dsdb = getdsdb;
fitdb = getfitdb;

% Offer to save session unless there's nothing to save
if isempty(down(dsdb)) && isempty(down(fitdb))
   resp = 'No';
else
  resp = questdlg('Save this Distribution Fitting session?', ...
                  'Distribution Fitting', 'Yes', 'No', 'Cancel', 'Yes');
end	
	 
if isempty(resp)
	resp = 'Cancel';
end

if isequal(resp,'Yes')
   ok = dfsession('save');
   if ~ok
      resp = 'Cancel';
   end
end

% Anything but cancel means go ahead and quit
if ~isequal(resp,'Cancel')
   set(gcbf,'CloseRequestFcn','');

   % Clear current session
   dfsession('clear');

   % Delete any dfittool-related figures
   h = dfgetset('evaluateFigure');
   if ~isempty(h) & ishandle(h), delete(h); end
   h = gcbf;
   if ~isempty(h) & ishandle(h), delete(h); end
end
