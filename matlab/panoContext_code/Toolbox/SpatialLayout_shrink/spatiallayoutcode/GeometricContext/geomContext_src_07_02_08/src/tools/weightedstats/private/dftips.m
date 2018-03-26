function dftips(varargin)
%DFTIPS Display data and fit tips for Distribution Fitting figure

%No input arguments are used here.

%   $Revision: 1.1.6.4 $  $Date: 2004/01/24 09:35:52 $
%   Copyright 2001-2004 The MathWorks, Inc.


dffig = gcbf;
if ~isequal(get(dffig,'SelectionType'),'normal')
   return;
end

h = hittest;
dftip = findobj(dffig,'tag','dftip');
dfdot = findobj(dffig,'tag','dfdot');
ftype = dfgetset('ftype');

% Figure out if the cursor is on something we know how to label
msg = '';
if (~isempty(h)) && ishandle(h) && isequal(get(h,'type'),'line')
   ax = get(h,'Parent');
   if isempty(ax) || ~isequal(get(ax,'Type'),'axes')
      ax = get(dffig,'CurrentAxes');
   end
   pt = get(ax,'CurrentPoint');
   x = pt(1,1);
   y = pt(1,2);
   htag = get(h,'tag');
   xlim = get(ax,'XLim');
   ylim = get(ax,'YLim');
   x = max(xlim(1), min(xlim(2),x));
   y = max(ylim(1), min(ylim(2),y));
   dx = diff(xlim) * 0.02;
   dy = 0;
   
   % Create a label for the the selected data set
   switch(htag)
    case {'dfdata' 'dfdbounds'}
      x = rounder(x,ax);  % make x a rounder number

      ds = get(h,'UserData');
      h = ds.line;
      xd = get(h,'XData');
      yd = get(h,'YData');
      if isequal(ftype,'probplot')
         % Label the closest point on the probability plot
         xyd = abs((xd-x)/diff(xlim)) + abs((yd-y)/diff(ylim));
         [ignore,jlo] = min(xyd);
         if ~isempty(jlo)
            jlo = jlo(1);
            x = xd(jlo);
            y = yd(jlo);
            cdffunc = getappdata(ax,'CdfFunction');
            plotparams = getappdata(ax,'DistributionParameters');
            yname = feval(cdffunc,y,plotparams{:});
            msg = sprintf('%s\n(%g, %g)',ds.name,x,yname);
         end

      elseif isequal(ftype,'pdf')
         % Label a bar of the histogram
         jlo = max(find(xd<=x));
         if ~isempty(jlo) && jlo<length(xd)
            xlo = xd(jlo);
            xhi = xd(jlo+1);
            y = yd(jlo);
            msg = sprintf(...
                '%s\nInterval [%g, %g]\nHeight = %g\nProbability = %g',...
                          ds.name, xlo, xhi, y, y*(xhi-xlo));
         end

      else
         % Label a cdf or similar empirical curve
         jlo = max(find(xd<=x));
         y = yd(jlo);
         dobounds = ds.showbounds;
         if dobounds
            if isempty(ds.boundline) || ~ishandle(ds.boundline)
               dobounds = false;
            end
         end
         if dobounds
            ybnds = get(ds.boundline,'YData');
            if isempty(ybnds);
               dobounds = false;
            else
               ylo = ybnds(jlo);
               yup = ybnds(length(yd) + 1 + jlo);
            end
         end
         if isempty(y)
            msg = '';
         elseif dobounds
            msg = sprintf('%s\nF(%g) = %g\nBounds [%g, %g]',...
                          ds.name,x,y,ylo,yup);
         else
            msg = sprintf('%s\nF(%g) = %g',ds.name,x,y);
         end
      end

    case {'distfit' 'dffbounds'}
      fit = get(h,'UserData');

      % Round x a bit
      if fit.iscontinuous || ~isequal(fit.ftype,'pdf')
         x = rounder(x,ax);
      else
         x = round(x);
      end

      descr = getdescription(fit);
      if isequal(ftype,'probplot')
         yname = eval(fit,x);
         invcdffunc = getappdata(ax,'InverseCdfFunction');
         plotparams = getappdata(ax,'DistributionParameters');
         y = feval(invcdffunc,yname,plotparams{:});
         msg = sprintf('%s\n%s\nF(%g) = %g)',fit.name,descr,x,yname);
      elseif ~isempty(fit.ybounds)
         [y,ylo,yup] = eval(fit,x);
         msg = sprintf('%s\n%s\nF(%g) = %g\nBounds [%g, %g]',...
                       fit.name,descr,x,y,ylo,yup);
      else
         y = eval(fit,x);
         msg = sprintf('%s\n%s\nF(%g) = %g',fit.name,descr,x,y);
      end
      dy = eval(fit,x+dx) - y;
   end
end

% If we can't label this thing, delete the label components
if isempty(msg)
   removetips(dffig);

% Otherwise we need to create the proper label
else
   if ~isempty(dfdot) && ishandle(dfdot)
      if isequal(x,get(dfdot,'XData')) && isequal(y,get(dfdot,'YData'))
         return;
      end
   end

   % Create the text and line components of the label if missing
   if isempty(dftip) || ~ishandle(dftip)
      yellow = [1 1 .85];
      dftip = text(x,y,'','Color','k','VerticalAlignment','bottom',...
                          'Parent',ax, 'Tag','dftip','Interpreter','none',...
                          'HitTest','off','FontWeight','bold',...
                          'Backgroundcolor',yellow,'Margin',3,...
                          'Edgecolor','k');
   end
   if isempty(dfdot) || ~ishandle(dfdot)
      dfdot = line(x,y,'Marker','o','LineStyle','none','Color','k',...
                       'Tag','dfdot','Parent',ax,'HitTest','off');
   end
   % Position the text so it is not clipped, then write the label
   if (x<sum(xlim)/2)
      ha = 'left';
   else
      ha = 'right';
      dx = - dx;
   end
   if isequal(htag,'dfdata')
      if y <= sum(ylim)/2
         va = 'bottom';
      else
         va = 'top';
      end
   else
      if sign(dx)==sign(dy)
         va = 'top';
      else
         va = 'bottom';
      end
   end

   set(dftip,'Position',[x+dx y 0],'String',msg,...
             'HorizontalAlignment',ha,'VerticalAlignment',va);
   set(dfdot,'XData',x,'YData',y);
   set(dffig, 'WindowButtonMotionFcn',@dftips);
   set(dffig, 'WindowButtonUpFcn',@disabletips);
end

% ------------------- remove tips from plot
function disabletips(varargin)
dffig = gcbf;
removetips(dffig);
set(dffig, 'WindowButtonUpFcn','');
set(dffig, 'WindowButtonMotionFcn','');

function removetips(dffig)
delete(findobj(dffig,'tag','dftip'));
delete(findobj(dffig,'tag','dfdot'));

% -------------------- round x before displaying tip there
function x = rounder(x,ax)

xlim = get(ax,'XLim');
xrange = diff(xlim);
if xrange>0
   if isequal(get(ax,'XScale'),'linear')
      pwr = floor(log10(0.005 * xrange));
      mult = 10 ^ -pwr;
      x = round(mult * x) / mult;
   elseif x>0
      pwr = floor(log10(x)) - 2;
      mult = 10 ^ -pwr;
      x = round(mult * x) / mult;
   end
end
