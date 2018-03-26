function dfupdatelegend(dffig,reset)
%DFUPLDATELEGEND Update legend in dfittool window

%   $Revision: 1.1.6.5 $  $Date: 2004/01/24 09:36:01 $ 
%   Copyright 2001-2004 The MathWorks, Inc. 

if nargin<2, reset=false; end

% If figure not passed in, find figure that contains this thing
while ~isequal(get(dffig,'parent'),0),
  dffig = get(dffig,'parent');
end

% Remember info about old legend, if any
ax = get(dffig,'CurrentAxes');
if isempty(get(ax,'Children'))
   legh = [];
else
   legh = legend('-find',ax);
end

% Try to put the legend where the curves are not likely to be.
% Position "0" is supposed to do this, but it is deprecated
% and can be slow.  Instead use a heurisitic.
ftype = dfgetset('ftype');
if ischar(ftype) && ~isempty(ftype) && ismember(ftype, {'pdf' 'survivor'})
   % The survivor function is decreasing.  Most pdf functions drop more
   % quickly to the right.  "top right" is probably good for these functions.
   legendpos = 'TR';
else
   % The remaining ones are increasing, so "top left" is probably good.
   legendpos = 'TL';
end
legend(ax, 'off');

% Maybe no legend has been requested
if isequal(dfgetset('showlegend'),'off')
   return
end

% Get data line handles and labels
hh = flipud(findobj(ax,'Type','line'));
hData = findobj(hh,'flat','Tag','dfdata');
n = length(hData);

textData = cell(n,1);
for j=1:length(hData)
   nm = '';
   ds = get(hData(j),'UserData');
   if ~isempty(ds) && ishandle(ds) && ~isempty(findprop(ds,'name'))
      nm = ds.name;
   end
   if isempty(nm)
      hData(j) = NaN;
   else
      textData{j} = nm;
   end
end
t = ~isnan(hData);
textData = textData(t);
hData = hData(t);
sortData = 1000*(1:length(hData));
if isempty(sortData)
   maxnum = 0;
else
   maxnum = max(sortData) + 1000;
end

% Indent bounds if there are two or more data set lines
if n>1
   pre = '  ';
else
   pre = '';
end

% Deal with confidence bounds, if any, around empirical cdf
n = length(hData);
textDataBounds = cell(n,1);
hDataBounds = repmat(NaN,n,1);
sortDataBounds = zeros(n,1);
for j=1:n
   ds = get(hData(j),'UserData');
   if ds.showbounds
      hbounds = ds.boundline;
      if ~isempty(hbounds) && ishandle(hbounds) ...
                           && ~isempty(get(hbounds,'YData'))
         textDataBounds{j} = [pre 'confidence bounds'];
         hDataBounds(j) = hbounds;
         sortDataBounds(j) = sortData(j) + .5;
      end
   end
end
if any(isnan(hDataBounds))
   t = isnan(hDataBounds);
   textDataBounds(t) = [];
   hDataBounds(t) = [];
   sortDataBounds(t) = [];
end

% Indent fits if there are two or more data lines
if (length(hData)>1)
   pre = '  ';
else
   pre = '';
end

% Get fit line handles and labels
hFit = findobj(hh,'flat','Tag','distfit');
sortFit = NaN*hFit;
n = length(hFit);
hFitConf = NaN*zeros(n,1);
textFit = cell(n,1);
nms = cell(n,1);
for j=1:length(hFit)
   try
      fit = get(hFit(j),'UserData');
      nm = fit.name;
   catch
      nm = '';
   end
   if isempty(nm)
      hFit(j) = NaN;
   else
      nms{j} = nm;
      textFit{j} = [pre nm];
      
      % Find the dataset for this fit
      ds = fit.dshandle;
      sortFitj = maxnum + j;
      for k=1:length(hData)
         if isequal(ds.name,textData{k})
            sortFitj = sortData(k) + j;
            break;
         end
      end
      sortFit(j) = sortFitj;

      % Look for bounds
      b = get(fit,'boundline');
      if ~isempty(b)
         hFitConf(j) = b(1);
      end
   end
end
t = ~isnan(hFit);
nms = nms(t);
textFit = textFit(t);
hFit = hFit(t);
sortFit = sortFit(t);
hFitConf = hFitConf(t);


% Indent bounds if there are two or more fits
if (length(hFit)>1)
   pre = [pre '  '];
end

% Get confidence bound line handles and labels
n = length(hFitConf);
textFitBounds = cell(n,1);
sortFitBounds = zeros(size(hFitConf));
for j=1:length(hFitConf)
   if ~isnan(hFitConf(j)) && ishandle(hFitConf(j)) ...
                          && ~isempty(get(hFitConf(j),'XData'))
      textFitBounds{j} = sprintf('%sconfidence bounds (%s)',pre,nms{j});
      sortFitBounds(j) = sortFit(j) + 0.5;
   else
      hFitConf(j) = NaN;
   end
end
t = ~isnan(hFitConf);
textFitBounds = textFitBounds(t);
hFitConf = hFitConf(t);
sortFitBounds = sortFitBounds(t);

% Combine everything together for the legend
h = [hData(:); hDataBounds(:); hFit(:); hFitConf(:)];
c = [textData; textDataBounds; textFit; textFitBounds];
s = [sortData(:); sortDataBounds(:); sortFit(:); sortFitBounds(:)];

% Sort so related things are together
[s,j] = sort(s);
c = c(j);
h = h(j);

% Create the legend
if (length(h)>0)
   ws = warning;
   lw = lastwarn;
   warning('off');
   oldu = get(ax,'Units');
   set(ax,'Units','normalized');
   try
      havecoords = isnumeric(legendpos) && length(legendpos)>1;
      if havecoords
          temppos = 0;
      else
          temppos = legendpos;
      end
      [legh,objh] = legend(ax,h,c,temppos);
      if havecoords
         set(legh,'Units','points');
         newpos = get(legh,'Position');
         newpos(1:2) = legendpos(1:2);
         set(legh,'Position',newpos);
      end
   catch
      objh = [];
   end
%   set(ax,'Units',oldu);
   warning(ws);
   lastwarn(lw);
   
   % Avoid treating ds/fit names as TeX strings
   objh = findobj(objh,'flat','Type','text');
   set(objh,'Interpreter','none','Hittest','off');
end

% Set a resize function that will handle legend and layout
set(dffig,'resizefcn','dfittool(''adjustlayout'');');
