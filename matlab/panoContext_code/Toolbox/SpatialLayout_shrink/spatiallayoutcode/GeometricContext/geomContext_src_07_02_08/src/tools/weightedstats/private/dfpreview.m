function [err, imsource] = dfpreview(dexpr, cexpr, fexpr, width, height, ds, binInfo)
% For use by DFITTOOL

%   $Revision: 1.1.6.6 $  $Date: 2004/01/24 09:35:44 $
%   Copyright 2001-2004 The MathWorks, Inc.

if nargin<6
   err = dfcheckselections(dexpr, cexpr, fexpr);
else
   err = '';
end
if ~isequal(err, '')
    imsource = [];
    return;
end

NONE='(none)';

if nargin<6
    data=evalin('base',dexpr);
else
    data = ds.y;
end

if nargin<6
    if isempty(cexpr) || isequal(cexpr, NONE)
        censoring=[];
    else
        censoring=evalin('base',cexpr);
    end
else
    censoring = ds.censored;
end

if nargin<6
    if isempty(fexpr) || isequal(fexpr, NONE)
        frequency=[];
    else
        frequency=evalin('base',fexpr);
    end
else 
    frequency=ds.frequency;
end

if nargin < 5
    width = 200;
    height = 200;
end

tempfigure=figure('units','pixels','position',[0 0 width height], ...
    'handlevisibility','callback', ...
    'integerhandle','off', ...
    'visible','off', ...
    'paperpositionmode', 'auto', ...
    'color','w');
tempaxes=axes('position',[.05 .05 .9 .9], ...
   'parent',tempfigure, ...
   'box','on', ...
   'visible','off');

% If data has a complex part, it will spit a warning to the command line, so
% turn off warnings before plotting.
warnstate=warning('off', 'all');


if nargin < 6
    binInfo = dfgetset('binDlgInfo');
elseif nargin < 7 
    binInfo = ds.binDlgInfo;
else
    % binInfo passed in
end

% If we're working on expressions rather than data in an existing data set,
% we may need to remove NaNs
[ignore1,ignore2,data,censoring,frequency] = statremovenan(data,censoring,frequency);

% Compute the bin centers using the ecdf
% to allow a quartile computation even when there is censoring.
[fstep, xstep] = ecdf(data, 'censoring', censoring, 'frequency', frequency);
[dum,binEdges] = dfhistbins(data,censoring,frequency,binInfo,fstep,xstep);

set(0,'CurrentFigure', tempfigure);
set(tempfigure,'CurrentAxes', tempaxes);

% Plot a histogram from ecdf using the computed number of bins
ecdfhist(tempaxes, fstep, xstep, 'edges', binEdges);
set(tempaxes, 'xtick',[],'ytick',[]);
axis(tempaxes,'tight');
allchildren = get(tempaxes, 'children');
patchchildren = findobj(allchildren,'flat','Type','patch');
set(patchchildren, 'facecolor', [.9 .9 .9]);
warning(warnstate);

x=hardcopy(tempfigure,'-dzbuffer','-r0');
% give the image a black edge
x(1,:,:)=0; x(end,:,:)=0; x(:,1,:)=0; x(:,end,:)=0;
imsource=im2mis(x);

delete(tempfigure);
