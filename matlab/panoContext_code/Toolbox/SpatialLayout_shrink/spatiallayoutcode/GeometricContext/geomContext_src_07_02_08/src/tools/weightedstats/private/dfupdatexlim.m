function dfupdatexlim(newminmax,updateplots)
%DFUPDATEXLIM Update the stored x axis min/max values

%   $Revision: 1.1.6.5 $  $Date: 2004/01/24 09:36:03 $
%   Copyright 2003-2004 The MathWorks, Inc.

minmax = [];                     % to become new x limits
oldminmax = dfgetset('xminmax'); % previous limits
ftype = dfgetset('ftype');
if nargin==0
   newminmax = [];
end

if isempty(newminmax) && isequal(ftype, 'icdf')
   % Default limits span most of the probability range
   minmax = [.01 .99];
elseif isempty(newminmax)
   % Update limits from datasets with a plotting flag on
   dsdb = getdsdb;
   ds = down(dsdb);
   while(~isempty(ds))
      if ds.plot == 1
         minmax = combineminmax(minmax,ds.xlim);
      end
      ds = right(ds);
   end

   % Update from fits with a plotting flag on
   fitdb = getfitdb;
   ft = down(fitdb);
   while(~isempty(ft))
      if ft.plot == 1
         minmax = combineminmax(minmax,xlim(ft));
      end
      ft = right(ft);
   end
else
   minmax = newminmax;
end

% Now update plot
dffig = dfgetset('dffig');
if ~isempty(minmax) && isequal(zoom(dffig,'getmode'),'off')
   ax = get(dffig,'CurrentAxes');
   islinscale = isequal(get(ax,'XScale'),'linear');
   if ~islinscale && any(minmax<=0)
      warning('stats:dfupdatexlim:NegativeDataIgnored',...
              'Negative data ignored.');
      minmax = [1e-6 1] * max(abs(minmax));
   end
   if isempty(newminmax) && ~isequal(ftype, 'icdf')
      % Adjust axis limits to include a margin around plotted points
      if islinscale
         dx = diff(minmax) * 0.01 * [-1 1];
         if all(dx==0), dx = [-1 1]; end
      else
         dlogx = .01 * diff(log(minmax));
         if dlogx==0, dlogx = 1; end
         dx = [minmax(1) * exp(-dlogx),  minmax(2) * exp(dlogx)] - minmax;
      end
   elseif minmax(1)==minmax(2)
      if islinscale
         dx = [-1 1];
      else
         dx = [minmax(1)/2, 2*minmax(1)];
      end
   else
      % Don't adjust the limits that were passed in or computed
      dx = 0;
   end
   oldxlim = get(ax,'XLim');
   newxlim = minmax + dx;
   if ~isequal(oldxlim,newxlim)
      set(ax,'XLim',newxlim);
      if nargin<2 || updateplots
         dfupdateallplots(false,true);
      end
   end
end
dfgetset('xminmax',minmax);

% ------------ Helper to combine old and new minmax values
function bothmm = combineminmax(oldmm,newmm)

if isempty(oldmm)
   bothmm = newmm;
elseif isempty(newmm)
   bothmm = oldmm;
else
   bothmm = [min(oldmm(1),newmm(1)) max(oldmm(2),newmm(2))];
end
