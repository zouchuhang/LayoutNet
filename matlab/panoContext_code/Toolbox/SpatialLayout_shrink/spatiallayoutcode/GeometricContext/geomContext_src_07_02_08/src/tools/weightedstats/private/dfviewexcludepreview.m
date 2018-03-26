function [imsource, D, C, F] = dfviewexcludepreview(outlier, width, height, dsname)
% For use by DFITTOOL

%   $Revision: 1.1.6.5 $
%   Copyright 2003-2004 The MathWorks, Inc.

tempfigure=figure('units','pixels','position',[0 0 width height], ...
      'handlevisibility','off', ...
      'integerhandle','off', ...
      'visible','off', ...
      'paperpositionmode', 'auto', ...
      'color','w');

NONE='(none)';

% We're excluding based on data in one dataset
ds = find(getdsdb,'name',dsname);
if isempty(ds)
   return
end

% Get data w/o NaNs but with no exclusion rule applied
[ydata, cens, freq] = getincludeddata(ds,[]);
ydata = real(ydata);
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

o = handle(outlier);
ylo = o.YLow;
if isempty(ylo)
	ylo = -Inf;
else
	ylo = str2double(ylo);
end

yhi = o.YHigh;
if isempty(yhi)
	yhi = Inf;
else
	yhi = str2double(yhi);
end

ylotest = o.YLowLessEqual;
yhitest = o.YLowLessEqual;

xlim = [min(x) max(x)];
xlim = xlim + .05 * [-1 1] * diff(xlim);
ylim = [min(y) max(y)];
ylim = ylim + .05 * [-1 1] * diff(ylim);
if ylim(1) == ylim(2)
    ylim = [0 2];
end    
ax=axes('position',[.05 .05 .9 .9], ...
      'parent',tempfigure, ...
      'xtick',[],'ytick',[], ...
      'box','on', ...
      'visible','off', 'XLim',xlim,'YLim',ylim);

if ylotest==0
   inbounds = x>=ylo;
else
   inbounds = x>ylo;
end
if yhitest==0
   inbounds = inbounds & x<=yhi;
else
   inbounds = inbounds & x<yhi;
end

figcolor = get(0,'defaultuicontrolbackgroundcolor');

t = inbounds;
l1 = line('XData',x(t),'YData',y(t),...
          'Color','b','Marker','.','LineStyle','none',...
          'Parent',ax);
t = ~inbounds;
l2 = line('XData',x(t),'YData',y(t),...
          'Color',figcolor/2,'Marker','.','LineStyle','none',...
          'Parent',ax);
alllines = [l1 l2];

gr = [.9 .9 .9];

% Create patches to show the area outside the domain and range

allpatches = [];

if (ylo ~= -Inf)
    xlo = max(ylo,xlim(1));
    xp = [xlim(1) xlo xlo xlim(1)];
    yp = [ylim(1) ylim(1) ylim(2) ylim(2)];
    p1=patch(xp,yp,gr,'LineStyle','none','Parent',ax);
    allpatches(end+1)= p1;
end

if (yhi ~= Inf)
    xhi = min(yhi,xlim(2));
    xp = [xlim(2) xhi xhi xlim(2)];
    yp = [ylim(1) ylim(1) ylim(2) ylim(2)];
    p2=patch(xp,yp,gr,'LineStyle','none','Parent',ax);
    allpatches(end+1)= p2;
end

set(ax,'Child',[alllines allpatches]);

x=hardcopy(tempfigure,'-dzbuffer','-r0');
% give the image a black edge
x(1,:,:)=0; x(end,:,:)=0; x(:,1,:)=0; x(:,end,:)=0;
imsource=im2mis(x);
delete(tempfigure);
[D, C, F] = dfviewdata(ds);