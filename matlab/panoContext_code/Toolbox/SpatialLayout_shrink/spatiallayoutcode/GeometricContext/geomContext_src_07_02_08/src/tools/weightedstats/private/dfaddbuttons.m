function dfaddbuttons(dffig)
%DFADDBUTTONS Add buttons to the curve fitting plot figure window

%   $Revision: 1.1.6.8 $  $Date: 2004/01/24 09:35:09 $
%   Copyright 2003-2004 The MathWorks, Inc.

% Clear out any old stuff
h0 = findall(dffig,'Type','uicontrol','Style','pushbutton');
if ~isempty(h0), delete(h0); end
p0 = ones(1,4);   % temporary position before adjustment

h0=uicontrol(dffig,'units','pixels','Tag','selectionframe',...
             'Style','frame','Position',p0);
h1=uicontrol(dffig,'units','pixels','Tag','displaytext',...
            'String','Display type:', 'style','text','Position',p0,...
            'HorizontalAlignment','left','FontWeight','bold');

choices = ['Density (PDF)|Cumulative probability (CDF)|Quantile (inverse CDF)|Probability plot|'...
           'Survivor function|Cumulative hazard'];
h2=uicontrol(dffig,'units','pixels','Tag','displaylist',...
            'String',choices, 'Style','pop','BackgroundColor',ones(1,3), ...
            'Callback', @cbkfunction,'Position',p0);
setappdata(h2,'codenames',...
             {'pdf' 'cdf' 'icdf' 'probplot' 'survivor' 'cumhazard'});

h3=uicontrol(dffig,'units','pixels','Tag','typetext',...
            'String','Distribution:', 'style','text','Position',p0,...
            'HorizontalAlignment','left','FontWeight','bold', 'Enable', 'off');

% Figure out which distributions could be used for probability plots,
% and remember their names and properties for later
alldist = dfgetdistributions;
dnames = {alldist.name};
dcodes = {alldist.code};
islocscale = [alldist.islocscale];
choices = dnames(islocscale);
default = strmatch('Normal',choices);
if length(default)~=1
   default = 1;
end
h4=uicontrol(dffig,'units','pixels','Tag','typelist',...
            'String',choices, 'Style','pop','BackgroundColor',ones(1,3),...
            'Value',default,'Enable','off', 'Callback',@cbkfunction,...
             'Position',p0);

% Store information about all distributions that could be used to make
% a probability plot, then separately list distributions that are ok
% given the current data.  Right now they're the same, but if we read
% in negative data we'll prune the "ok" list to omit positive
% distributions.
setappdata(h4,'allfullnames',choices);
setappdata(h4,'allcodenames',dcodes(islocscale));
setappdata(h4,'alldistspec',alldist(islocscale));
setappdata(h4,'okcodenames',dcodes(islocscale));

hvec = [h0 h1 h2 h3 h4];
setappdata(dffig,'selectioncontrols',hvec);

% Define information for the buttons
strings = {xlate('Data...') xlate('New Fit...') xlate('Manage Fits...') xlate('Evaluate...') xlate('Exclude...')};
tips = {xlate('Import, view, rename, plot and delete data') ...
        xlate('Add a fitted distribution') ...
        xlate('Edit, view, plot and rename fits') ...
        xlate('Evaluate fits to compute a table of results') ...
        xlate('Define rules for excluding data from a fit')};
cbacks = {@cbkdata @cbknewfit @cbkmanfit @cbkeval @cbkexclude};

% Add the buttons to the figure
tags = {'dfdata' 'dfnewfit' 'dfmanfit' 'dfevaluate' 'dfexclude'};
n = length(strings);
h = zeros(1,n);
for j=1:length(strings)
   h(j) = uicontrol(dffig,'Units','pixels', ...
                    'Position',[j,1,1,1],...
                    'String',strings{j}, 'TooltipString',tips{j}, ...
                    'Callback',cbacks{j}, 'Tag',tags{j});
end
setappdata(dffig,'buttoncontrols',h);

% ---------------------- callback for Import button
function cbkdata(varargin)
%CBKDATA Callback for Data button

delete(findall(gcbf,'Tag','dfstarthint'));
com.mathworks.toolbox.stats.Data.showData;

% ---------------------- callback for New Fit button
function cbknewfit(varargin)
%CBKNEWFIT Callback for New Fit button

com.mathworks.toolbox.stats.Fitting.getFitting.showNewFit;


% ---------------------- callback for Manage Fit button
function cbkmanfit(varargin)
%CBKMANFIT Callback for Manage Fit button

com.mathworks.toolbox.stats.FitManager.showFitManager;

% ---------------------- callback for Evaluate button
function cbkeval(varargin)
%CBKEVAL Callback for Evaluate button

% Get the current horizontal axis limits of the main figure.  We'll set the
% default points at which to evaluate to a "pretty" colon expression that
% spans those limits, with 10 steps.
xlims = dfgetset('xminmax');

if isempty(xlims)
   xstr = '';
else
   % Start out by choosing a rounding that will give the smaller (in
   % magnitude) of min(x) and max(x) a single sig digit, but at most two sig
   % digits in the larger.  The smaller may round to zero.  If the min and max
   % round to the same thing, use more digits until we get rounded numbers
   % that differ.
   xmag = max(min(abs(xlims)), max(abs(xlims))/10);
   rounder = 10^floor(log10(xmag));
   while true
       xmin = round(xlims(1)./rounder) * rounder;
       xmax = round(xlims(2)./rounder) * rounder;
       if xmin < xmax, break; end
       rounder = rounder/10;
   end
   
   % Create 10 steps, where the step size will have three more significant
   % digits than the endpoints.
   stepRounder = rounder ./ 1000;
   step = floor((xmax-xmin)./(10*stepRounder)) * stepRounder;
   
   % Figure out how many digits we need display in order to distinguish the
   % endpoints.  That's the number of digits to the left of the decimal, plus
   % however many we've rounded to on the right.  Allow for at least four so
   % that we'll get, e.g., "4000", and not "4e+03".
   xminDigits = max(round(log10(max(abs(xmin),1))-log10(rounder)),4);
   xmaxDigits = max(round(log10(max(abs(xmax),1))-log10(rounder)),4);
   xstr = sprintf('%0.*g:%g:%0.*g', xminDigits, xmin, step, xmaxDigits, xmax);
end

% Get the current plot type of the main figure, we'll set the default
% function type to evaluate based on that.
ftype = dfgetset('ftype');
if strcmp(ftype,'probplot')
    ftype = 'cdf'; % can't evaluate a prob plot, use cdf instead
% else {'pdf' 'cdf' 'survivor' 'icdf' 'cumhazard' 'hazrate'}
    % otherwise use the current setting
end

com.mathworks.toolbox.stats.Evaluate.showEvaluate(ftype,xstr);

% ---------------------- callback for Exclude button
function cbkexclude(varargin)
%CBKEXCLUDE Callback for Exclude button

com.mathworks.toolbox.stats.Exclude.showExclude;

% ---------------------- callback for display list
function cbkfunction(varargin)
%CBKFUNCTION Callback for setting function to display

% Get the requested function and distribution types
dffig = gcbf;
dfsetplottype(dffig);
