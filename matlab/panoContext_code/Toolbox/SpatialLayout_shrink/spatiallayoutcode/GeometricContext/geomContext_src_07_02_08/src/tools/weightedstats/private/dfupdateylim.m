function dfupdateylim
%DFUPDATEYLIM Update the y axis min/max values

%   $Revision: 1.1.6.4 $  $Date: 2004/01/24 09:36:04 $
%   Copyright 2003-2004 The MathWorks, Inc.

dminmax = [];                    % to indicate y data limits

% Check y limits of all fits
fminmax = [];
fitdb = getfitdb;
ft = down(fitdb);
while(~isempty(ft))
   if ft.plot==1 && ~isempty(ft.linehandle) && ishandle(ft.linehandle)
      fminmax = combineminmax(fminmax,ft.ylim);
   else
      ds = ft.dshandle;
      if ds.plot==1 && ~isempty(ds.line) && ishandle(ds.line)
         fminmax = combineminmax(fminmax,ds.ylim);
      end
   end
   ft = right(ft);
end

% Check any datasets with a plotting flag on
dsdb = getdsdb;
ds = down(dsdb);
while(~isempty(ds))
   if ds.plot == 1
      dminmax = combineminmax(dminmax,ds.ylim);
   end
   ds = right(ds);
end

% Adjust data min/max to take fits into account, but don't allow
% fit extrapolations to overwhelm data values
if isempty(dminmax)
   if isempty(fminmax), return; end
   dminmax = fminmax;
elseif ~isempty(fminmax)
   dy = diff(dminmax);
   dminmax(1) = max(dminmax(1)-dy/2, min(dminmax(1),fminmax(1)));
   dminmax(2) = min(dminmax(2)+dy/2, max(dminmax(2),fminmax(2)));
end

dffig = dfgetset('dffig');
ax = get(dffig,'CurrentAxes');
if isequal(get(ax,'YScale'),'linear')
   dy = diff(dminmax) * 0.05 * [-1 1];
   if all(dy==0)
      dy = [-1 1];
   end
elseif dminmax(2)>dminmax(1)
   dlogy = .01 * diff(log(dminmax));
   if (dlogy==0), dlogy = 1; end
   dy = [dminmax(1) * exp(-dlogy),  dminmax(2) * exp(dlogy)] - dminmax;
else
   dy = 0;
end
   
ftype = dfgetset('ftype');
if isempty(ftype)   % may happen during initialization
   ftype = 'pdf';
end
switch(ftype)
 case {'cdf' 'survivor'}
   % Bounded functions, no need to extend in either direction
   dy(:) = 0;
   
 case {'pdf' 'cumhazard'}
   % Positive functions, no need to extend below zero
   dy = max(0,dy);
end

if isequal(zoom(dffig,'getmode'),'off')
   set(ax,'YLim',dminmax+dy);
end


% ------------ Helper to combine old and new minmax values
function bothmm = combineminmax(oldmm,newmm)

if isempty(oldmm)
   bothmm = newmm;
elseif isempty(newmm)
   bothmm = oldmm;
else
   bothmm = [min(oldmm(1),newmm(1)) max(oldmm(2),newmm(2))];
end
