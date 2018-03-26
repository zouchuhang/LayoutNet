function dfupdateppdists(dffig)
%DFUPDATEPPDISTS Update distribution list for probability plots

%   $Revision: 1.1.6.5 $  $Date: 2004/01/24 09:36:02 $
%   Copyright 2003-2004 The MathWorks, Inc.

if nargin<1 || isempty(dffig)
   dffig = dfgetset('dffig');
end

% Get handle to control containing the distribution list
hsel = getappdata(dffig,'selectioncontrols');

% Determine plot type, don't continue unless it's a probability plot
h = hsel(3);   % handle of display type control
choice = get(hsel(3),'Value');
ftypes = getappdata(hsel(3),'codenames');
ftype = ftypes{choice};
if ~isequal(ftype, 'probplot');
   return
end

% Look at all plotted data sets and check for negative data
dsdb = getdsdb;
dslist = find(dsdb);
dslist(dslist==dsdb) = [];
nds = length(dslist);
xmin = 1;
for j=1:nds
   ds = dslist(j);
   if ds.plot
      xlimits = xlim(ds);
      if (xlimits(1) < xmin)
         xmin = xlimits(1);
         if (xmin<0)
            break
         end
      end
   end
end

% Get the current selection, and try to re-select it later
h = hsel(5);   % handle of distribution control
dlist = get(h, 'String');
choice = get(h, 'Value');
if choice<=length(dlist)
   cursel = dlist{choice};
else
   cursel = [];
end

% Get the distribution entries only
dlist = getappdata(h, 'allfullnames');
codenames = getappdata(h, 'allcodenames');
ndists = length(dlist);

% With 0 or negative data showing, omit some distributions
if (xmin<=0)
   ok = true(ndists,1);
   allspec = getappdata(h,'alldistspec');
   
   % ------- Could vectorize the following
   for j=1:ndists
      spec = allspec(j);
      ok(j) = checkdist(spec,xmin);
   end
   % -------
   
   dlist = dlist(ok);
   codenames = codenames(ok);
   ndists = length(dlist);
end
setappdata(h,'okcodenames',codenames);
   
% Create a list of fit names
fitdb = getfitdb;
flist = find(fitdb);
flist(flist==fitdb) = [];
nfits = length(flist);

% Weed out fits that cannot be used for probability plotting
for j=nfits:-1:1
   fj = flist(j);
   if ~fj.isgood || isequal(fj.fittype,'smooth')
      flist(j) = [];
   else
      distspec = fj.distspec;
      if isempty(distspec) || ~distspec.iscontinuous ...
                           || ~checkdist(distspec,xmin)
         flist(j) = [];
      end
   end
end
nfits = length(flist);

% Create combined list
newdlist = cell(ndists + nfits, 1);
newdlist(1:ndists) = dlist(:);
savelist = cell(nfits,1);
for j=1:nfits
   fj = flist(j);
   newdlist{ndists+j} = sprintf('Estimated %s (%s)',fj.distspec.name,fj.name);
   savelist{j} = fj.name;
end
setappdata(h, 'fitnames', savelist);

% Re-select the previous selection
set(h, 'String', newdlist);
choice = 1;
if ~isempty(cursel)
   choices = strmatch(cursel,newdlist,'exact');
   if numel(choices) == 1
      choice = choices;
   end
end
set(h, 'Value', choice);


% -----------------------------
function ok = checkdist(spec,xmin)
%Check distribution against minimum data value

lobnd = spec.support(1);
strict = ~spec.closedbound(1);
if strict && lobnd>=xmin
   ok = false;
elseif ~strict && lobnd>xmin
   ok = false;
else
   ok = true;
end
