function dfgraphexclude(dsname,xlo,xhi)
%DFGRAPHEXCLUDE  Create graph for selecting (x,y) pairs to exclude
%   DFGRAPHEXCLUDE(EXCLUDEPANEL,DSNAME,LOBND,UPBND) creates a graph
%   tied to the Java exclusion panel EXCLUDEPANEL, for dataset DSNAME, with
%   current lower and upper bounds LOBND and UPBND.  It provides a graphical
%   way to modify those bounds.

% Copyright 2001-2004 The MathWorks, Inc.
% $Revision: 1.1.6.8 $  $Date: 2004/02/01 22:10:38 $ 


% Use old figure if any, or create a new one with a plot of the data
t = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');
c = get(0,'Child');
f = findobj(c,'flat','Type','figure','Tag','dfexcludegraph');
set(0,'ShowHiddenHandles',t);
if ~isempty(f)
   subfig = f;
else
   subfig = setupfigure(dsname,xlo,xhi);
end
if isempty(subfig)
   return
end

figure(subfig)

% Adjust the patches to show the desired
ax = get(subfig,'CurrentAxes');
xlim = get(ax,'XLim');

% If bounds already exist, put them onto the graph
if nargin>=2 && ~isempty(xlo) && ~isinf(xlo) && xlo > xlim(1)
   addremovepatch(subfig,'lower',xlo,'add')
end
if nargin>=3 && ~isempty(xhi) && ~isinf(xhi) && xhi < xlim(2)
   addremovepatch(subfig,'upper',xhi,'add')
end

set(ax,'XLimMode','manual','YLimMode','manual');
set(subfig,'HandleVisibility','callback');
dfgetset('dfsubfig',subfig);

return


% -------------- Create figure
function subfig = makefigure(x,y,dolegend,xlo,xhi)

figcolor = get(0,'defaultuicontrolbackgroundcolor');
subfig = figure('IntegerHandle','off','Units','pixels',...
                'HandleVisibility','callback',...
                'name','Define Boundary for Exclusion Rule',...
                'numbertitle','off',...
                'color',figcolor,...
                'Tag','dfexcludegraph',...
                'DeleteFcn', @closefig,...
                'doublebuffer','on',...
                'Dock','off');

% Remove menus
set(subfig,'Menubar','none');

% Restore toolbar but keep only zoom tools
set(subfig,'toolbar','figure');
h = findall(subfig,'Type','uitoolbar');
h1 = findall(h);        % Get all children
h1(h1==h) = [];         % Not including the toolbar itself
h2 = findall(h1,'flat','TooltipString','Zoom In');
h1(h2==h1) = [];
h2 = findall(h1,'flat','TooltipString','Zoom Out');
h1(h2==h1) = [];
delete(h1);

% Add axes, for now in default position
xlim = [min(x) max(x)];
if isfinite(xlo)
    xlim(1) = min(xlim(1), xlo);
end
if isfinite(xhi)
    xlim(2) = max(xlim(2), xhi);
end

xlim = xlim + .1 * [-1 1] * diff(xlim);
ax = axes('Parent',subfig,'Box','on','HitTest','off',...
          'XLim',xlim,'YLim',[0, max(y)+1]);

% Place buttons where required
p = get(subfig,'Position');

pbutton = [5 5 50 20];
a1=uicontrol(subfig,'Units','pixels','Position',pbutton,...
            'Tag','close','Callback',@done,...
            'String','Close');
e = get(a1,'Extent');
pbutton(3:4) = 1.5 * e(3:4);
set(a1,'Position',pbutton);

margin = 15;

% Use longer string to get extent, then replace with shorter one
a3 = uicontrol(subfig,'Units','pixels','Position',pbutton,...
              'Tag','upper','Callback',@buttoncallback,...
              'String','Remove Upper Limit');
e = get(a3,'Extent');
pbutton(3) = 1.25 * e(3);
pbutton(1) = max(1,pbutton(1) - margin - pbutton(3));
set(a3,'Position',pbutton,'String','Add Upper Limit');

pbutton(1) = max(5, pbutton(1)-margin-pbutton(3));
a4=uicontrol(subfig,'Units','pixels','Position',pbutton,...
            'Tag','lower','Callback',@buttoncallback,...
            'String','Add Lower Limit');
setappdata(subfig,'buttons',[a1 a3 a4]);

% Place text as well
if dolegend
   p = [5 5 10 10];
   a = uicontrol(subfig, 'Units','pixels', 'Position',p,...
                 'Style','text','String','Observed = filled circles');
   e = get(a,'Extent');
   p(3:4) = e(3:4);
   set(a,'Position',p);
   p(2) = p(2) + p(4) + 5;
   a = uicontrol(subfig, 'Units','pixels', 'Position',p,...
                 'Style','text','String','Censored = open circles');
end

% Leave axis units as pixels, and finally set resize function
set(ax,'Units','pixels');
resize(subfig);
set(subfig,'ResizeFcn',@resize,'WindowButtonMotionFcn',@fixcursor);

% ------------------- helper function to set up figure
function subfig = setupfigure(dsname,xlo,xhi)
%SETUPFIGURE Set up figure to do graphical exclusion
% We're excluding based on data in one dataset
dsdb = dfswitchyard('getdsdb');
a = down(dsdb);
ds = [];
while(~isempty(a))
   if isequal(dsname,a.name)
      ds = a;
      break;
   end
   a = right(a);
end
if isempty(ds)
   subfig = [];
   return
end
   
[ydata,cens,freq] = getincludeddata(ds,[]); % get data w/o NaNs
if isempty(cens)
   cens = zeros(size(ydata));
end
if isempty(freq)
   freq = ones(size(ydata));
end

% Sort y and carry along the rest
[ydata,i] = sort(ydata);
cens = cens(i);
freq = freq(i);

% Create x and y vectors to plot
n = sum(freq);
x = zeros(n,1);
y = zeros(n,1);
g = zeros(n,1);
j = 1;
x(1:freq(1)) = ydata(1);
y(1:freq(1)) = (1:freq(1))';
g(1:freq(1)) = cens(1);
i = freq(1)+1;
for k=2:length(ydata)
   for j=1:freq(k)
      x(i) = ydata(k);
      g(i) = cens(k);
      if (i>1) && (x(i)==x(i-1))
         y(i) = y(i-1) + 1;
      else
         y(i) = 1;
      end
      i = i+1;
   end
end

% Make a figure to receive graph
dolegend = any(g==0) & any(g==1);
subfig = makefigure(x,y,dolegend, xlo, xhi);
ax = get(subfig,'CurrentAxes');

% Place data points into graph
t = (g==0);
if any(t)
   line('XData',x(t),'YData',y(t),'HitTest','off',...
          'Color','b','Marker','.','LineStyle','none',...
          'MarkerSize',24,'Parent',ax,'Tag','observed');
end
t = (g==1);
if any(t)
   line('XData',x(t),'YData',y(t),'HitTest','off',...
          'Color','b','Marker','o','LineStyle','none',...
          'Parent',ax,'Tag','censored');
end

% -----------------------
function resize(subfig,varargin)

if ~isempty(subfig) && ~ishandle(subfig) ...
                    && ~isequal(get(subfig,'type'),'figure')
   subfig = gcbf;
end
pFig = get(subfig,'Position');
ax = get(subfig,'CurrentAxes');
oldpos = getappdata(subfig,'oldpos');
if isequal(oldpos,pFig)
   return
end

% Position buttons against the edge
paxes = get(ax,'Position');
margin = 15;
hButtons = getappdata(subfig,'buttons');
base = pFig(3);
for j=1:length(hButtons)
   p = get(hButtons(j),'Position');
   p(1) = base - p(3) - margin;
   set(hButtons(j),'Position',p);
   base = p(1);
end

% Move axes above buttons
paxes(1) = 75;
paxes(3) = max(1,pFig(3)-150);
paxes(2) = p(2) + p(4) + 3*margin;
paxes(4) = max(1,pFig(4) - margin - paxes(2));
set(ax,'Position',paxes);
setappdata(subfig,'oldpos',pFig);

% --------------------------------
function buttoncallback(varargin)

% Get some handles and dimensions
button = gcbo;
fig = gcbf;
ax = get(fig,'CurrentAxes');
xlim = get(ax,'XLim');
ylim = get(ax,'YLim');
dx = 0.05 * diff(xlim);
buttontag = get(button,'Tag');
addremovepatch(fig,buttontag);
updateGUI;

% --------------------------------
function addremovepatch(fig,whichbound,xbnd,addremove)

% Get some handles and dimensions
ax = get(fig,'CurrentAxes');
xlim = get(ax,'XLim');
ylim = get(ax,'YLim');
dx = 0.05 * diff(xlim);
if nargin>=4 && isequal(addremove,'add')
   forceadd = true;
else
   forceadd = false;
end

if nargin<3
   if isequal(whichbound,'lower')
      xbnd = xlim(1) + dx;
   else
      xbnd = xlim(2) - dx;
   end
end

% Carry out requested action
if isequal(whichbound,'lower')
   hPatch = findall(fig,'Tag','lowerpatch');
   if isempty(hPatch) || forceadd
      otherpatch = findall(fig,'Tag','upperpatch');
      if ~isempty(otherpatch)
         % Never put new limit beyond the other limit
         otherx = get(otherpatch,'XData');
         otherx = otherx(2);
         xbnd = min(xbnd, xlim(1) + .9*(otherx-xlim(1)));
      end
      x = [xlim(1), xbnd,    xbnd,    xlim(1), xlim(1)];
      y = [ylim(1), ylim(1), ylim(2), ylim(2), ylim(1)];
      if isempty(hPatch)
         patch(x,y,[.9 .9 .9],'Parent',ax,'Tag','lowerpatch',...
            'FaceAlpha',0.6,'ButtonDownFcn',@startselect);
      else
         set(hPatch,'XData',x,'YData',y);
      end
         
      newtxt = 'Remove Lower Limit';
   else
      delete(hPatch);
      newtxt = 'Add Lower Limit';
   end
else
   hPatch = findall(fig,'Tag','upperpatch');
   if isempty(hPatch) || forceadd
      otherpatch = findall(fig,'Tag','lowerpatch');
      if ~isempty(otherpatch)
         % Never put new limit beyond the other limit
         otherx = get(otherpatch,'XData');
         otherx = otherx(2);
         xbnd = max(xbnd, xlim(2) - .9*(xlim(2)-otherx));
      end
      x = [xlim(2), xbnd,    xbnd,    xlim(2), xlim(2)];
      y = [ylim(1), ylim(1), ylim(2), ylim(2), ylim(1)];
      if isempty(hPatch)
         patch(x,y,[.9 .9 .9],'Parent',ax,'Tag','upperpatch',...
            'FaceAlpha',0.6,'ButtonDownFcn',@startselect);
      else
         set(hPatch,'XData',x,'YData',y);
      end
      newtxt = 'Remove Upper Limit';
   else
      delete(hPatch);
      newtxt = 'Add Upper Limit';
   end
end

% Update button text
button = findobj(fig,'Tag',whichbound);
set(button,'String',newtxt);

% ------------- function to initiate graphical selection
function startselect(varargin)

% Get figure and axis handles, define functions to do and end selection
subfig = gcbf;
ax = get(subfig,'CurrentAxes');

% Get current exclusion limits, use axis limits if none
lims = get(ax,'XLim');
hPatch = findall(subfig,'Tag','lowerpatch');
if ~isempty(hPatch)
   x = get(hPatch,'XData');
   lims(1) = max(x(:));
end
hPatch = findall(subfig,'Tag','upperpatch');
if ~isempty(hPatch)
   x = get(hPatch,'XData');
   lims(2) = min(x(:));
end

% Save information for other functions
hPatch = gcbo;
set(subfig,'WindowButtonMotionFcn',{@movepatch hPatch},...
           'WindowButtonUpFcn',@endmove);
setappdata(ax,'limits',lims);
setappdata(ax,'objmoving',hPatch);


% ------------- function to update GUI
function updateGUI()
subfig = gcbf;

hPatch = findall(subfig,'Tag','lowerpatch');
if isempty(hPatch)
   xl = '';
else
   xl = get(hPatch,'XData');
   xl = num2str(xl(2));
end
hPatch = findall(subfig,'Tag','upperpatch');
if isempty(hPatch)
   xh = '';
else
   xh = get(hPatch,'XData');
   xh = num2str(xh(2));
end

com.mathworks.toolbox.stats.Exclude.getExcludePanel.updateBoundsFields(xl, xh);


% ------------- function to complete graphical selection
function endmove(varargin)

% Turn off window functions to end selection
subfig = gcbf;
set(subfig,'WindowButtonMotionFcn',@fixcursor, 'WindowButtonUpFcn',[]);
updateGUI;

% ------------- move patch boundary
function movepatch(ignore1,ignore2,hPatch)

varargin
subfig = gcbf;
ax = get(subfig,'CurrentAxes');

% Get exclusion limits and axis limits
lims = getappdata(ax,'limits');
xlim = get(ax,'XLim');
delta = .01 * abs(xlim(2) - xlim(1));

% Extend patch to the current point, but within limits
cp = get(gca,'CurrentPoint');
x = cp(1);
if isequal(get(hPatch,'Tag'),'lowerpatch')
   lobnd = xlim(1) + min(delta, .9*(lims(2)-xlim(1)));
   upbnd = min(xlim(2),lims(2)) - delta;
   x = max(lobnd, min(x,upbnd));
   lims(1) = x;
else
   lobnd = max(xlim(1),lims(1)) + delta;
   upbnd = xlim(2) - min(delta, .9*(xlim(2)-lims(1)));
   x = min(upbnd, max(lobnd,x));
   lims(2) = x;
end

% Update saved limits, and x data for this patch
setappdata(ax,'limits',lims);
xdata = get(hPatch,'XData');
xdata(2:3) = x;
set(hPatch,'XData',xdata);

% --------------- set cursor if we're on something that can move
function fixcursor(varargin)
ptr = get(gcbf,'Pointer');
onpatch = isequal(get(hittest,'Type'),'patch');
if isequal(ptr,'arrow')
    if onpatch
        set(gcbf,'Pointer','left');
    end
else
    if ~onpatch
        set(gcbf,'Pointer','arrow');
    end
end

% --------------- close figure
function done(varargin)

delete(gcbf);

% ---------------------- helper to notify GUI that figure is closing
function closefig(varargin)
%CLOSEFIG 

com.mathworks.toolbox.stats.Exclude.getExcludePanel.setGraphExcludeFlag(false);
closereq;

