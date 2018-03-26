function dfsetfunction(dffig,ftype,dtype)
%DFSETFUNCTION Set value of probability function to be plotted

%   $Revision: 1.1.6.9 $  $Date: 2004/01/24 09:35:50 $
%   Copyright 2003-2004 The MathWorks, Inc.

% Figure out what we're changing from
oldftype = dfgetset('ftype');
olddtype = dfgetset('dtype');

% Remember these settings
dfgetset('ftype',ftype);
if nargin>=3
   dfgetset('dtype',dtype);
else
   dtype = [];
end

% Get rid of axis controls
if isequal(dfgetset('showaxlimctrl'),'on')
   dftoggleaxlimctrl(dffig);
end

% Determine and remember if this function type supports bounds
oktypes = {'cdf' 'survivor' 'cumhazard' 'icdf'};
dobounds = ismember(ftype, oktypes);
dfgetset('dobounds',dobounds);

% Get the array of data sets
dsdb = getdsdb;
ds = down(dsdb);

if ~isequal(oldftype,ftype)
   % Change the function type for each one
   while(~isempty(ds))
      setftype(ds,ftype);
      ds = right(ds);
   end
elseif isequal(ftype,'probplot') && ~isequal(olddtype,dtype)
   while(~isempty(ds))
      clearplot(ds);
      ds = right(ds);
   end
end   

% Get the array of fits
fitdb = getfitdb;
ft = down(fitdb);
referenceFit = [];

if ~isequal(oldftype,ftype)
   % Change the function type for each one
   while(~isempty(ft))
      setftype(ft,ftype);
      ft = right(ft);
   end
end

% Determine if a specific set of parameters (reference fit) is required
if isequal(ftype,'probplot')
   if ishandle(dtype)
      referenceFit = dtype;
   else
      referenceFit = [];
   end
end

% Update title, labels, appdata, and (for probability plots) axes
ax = get(dffig,'CurrentAxes');
newplot(ax);
setappdata(ax,'ReferenceDistribution','');
setappdata(ax,'CdfFunction','');
setappdata(ax,'InverseCdfFunction','');
setappdata(ax,'DistributionParameters','');
setappdata(ax,'LogScale','');

% Define the colors to be used here
a = [3 0 2 1 3 3 3 2 2 0 2 3 0 1 2 1 0 1 0 1 1
     0 0 1 1 0 3 2 2 1 2 0 1 3 2 3 0 1 3 0 2 0
     0 3 0 1 3 0 1 2 3 1 1 2 0 3 1 2 2 2 0 0 2]'/3;
set(ax,'ColorOrder',a);

% Turn on grid if requested
if isequal(dfgetset('showgrid'), 'on')
   set(ax,'xgrid','on','ygrid','on')
end

if isequal(ftype,'probplot')
   if isempty(referenceFit)
      probplot(ax,dtype);
   else
      probplot(ax,{referenceFit.distspec, referenceFit.params})
   end
   title(ax,'');
elseif isequal(ftype,'icdf')
   set(get(ax,'XLabel'),'String','Probability');
   set(get(ax,'YLabel'),'String','Quantile');
else
   othertypes =  {'pdf'               'cdf'                    ...
                  'survivor'          'cumhazard'};
   otherlabels = {'Density'           'Cumulative probability' ...
                  'Survivor function' 'Cumulative hazard'};
   jtype = strmatch(lower(ftype),othertypes,'exact');
   if isempty(jtype)   % should never happen
      ylab = ftype;
   else
      ylab = otherlabels{jtype};
   end
   set(get(ax,'XLabel'),'String','Data');
   set(get(ax,'YLabel'),'String',ylab);
end
set(ax, 'box','on','Tag','main',...
        'XLimMode','manual','YLimMode','manual','ZLimMode','manual',...
        'CLimMode','manual','AlimMode','manual');

% Reset the x limits, update plotted curves, and set y limits
dfupdateallplots(true,false);    % update data sets
dfupdatexlim([],false);          % get new x limits
dfupdateallplots(false,true);    % update fits
dfupdateylim;                    % now compute y limits

% Update the legend
dfupdatelegend(dffig);

% Make sure each data set's plot property is enabled as appropriate
if isequal(oldftype,'probplot') || isequal(ftype,'probplot')
   ds = down(dsdb);
   while(~isempty(ds))
      com.mathworks.toolbox.stats.DataSetsManager.getDataSetsManager.dataSetChanged(...
          java(ds), ds.name, ds.name);
      ds = right(ds);
   end
end