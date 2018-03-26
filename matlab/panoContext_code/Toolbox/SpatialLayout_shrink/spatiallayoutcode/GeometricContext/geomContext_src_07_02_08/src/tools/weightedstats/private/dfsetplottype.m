function dfsetplottype(dffig,ftype,dtype)
%DFSETPLOTTYPE Set plot type and distribution for distribution fitting

%   $Revision: 1.1.6.4 $  $Date: 2004/01/24 09:35:51 $
%   Copyright 2003-2004 The MathWorks, Inc.

hFunctionList = findall(dffig, 'Tag', 'displaylist');
hDistributionList = findall(dffig, 'Tag', 'typelist');

% Determine plot type
if nargin<2 || isempty(ftype)
   choice = get(hFunctionList,'Value');
   ftypes = getappdata(hFunctionList,'codenames');
   ftype = ftypes{choice};
end
isprobplot = isequal(ftype, 'probplot');

% Make sure distribution list is accurate
dfupdateppdists(dffig);

% Get distribution type for probability plot
if isprobplot && (nargin<3 || isempty(dtype))
   choice = get(hDistributionList,'Value');
   dtypes = getappdata(hDistributionList,'okcodenames');
   ntypes = length(dtypes);
   
   % Choose either a distribution name or a fit
   if choice<=ntypes
      dtype = dtypes{choice};
   else
      flist = getappdata(hDistributionList, 'fitnames');
      fname = flist{choice-ntypes};
      dtype = find(getfitdb, 'name', fname);
   end
elseif ~isprobplot
   dtype = [];
end
   
% Enable or disable the distribution field for probability plots
typetext = findall(dffig, 'tag', 'typetext');
if isprobplot
    set(hDistributionList, 'Enable', 'on');
    set(typetext, 'Enable', 'on');
else
    set(hDistributionList, 'Enable', 'off');
    set(typetext, 'Enable', 'off');
end

% Update everything to use this type
dfsetfunction(dffig,ftype,dtype);
