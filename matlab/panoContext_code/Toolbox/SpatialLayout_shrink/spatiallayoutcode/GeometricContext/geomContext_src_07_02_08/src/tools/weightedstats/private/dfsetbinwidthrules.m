function dfsetbinwidthrules(dataset, v0, v1, v2, v3, v4, v5, v6)
% DFSETBINWIDTHRULES Helper function for the dfittool set default bin width panel

%   Copyright 2003-2004 The MathWorks, Inc. 
%   $Revision: 1.1.6.3 $  $Date: 2004/01/24 09:35:47 $

binDlgInfo.rule = v0;
binDlgInfo.nbinsExpr = v1;
binDlgInfo.nbins = str2num(v1);
binDlgInfo.widthExpr = v2;
binDlgInfo.width = str2num(v2);
binDlgInfo.placementRule = v3;
binDlgInfo.anchorExpr = v4;
binDlgInfo.anchor = str2num(v4);
binDlgInfo.applyToAll = v5;
binDlgInfo.setDefault = v6;


if (v5 == true) % apply to all
    dsdb = getdsdb;
    ds = down(dsdb);
    while(~isempty(ds))
        ds.binDlgInfo = binDlgInfo;
        ds = right(ds);
    end
    dfupdateallplots(true, false, true); 
elseif ~isempty(dataset)
    ds = handle(dataset);
    ds.binDlgInfo = binDlgInfo;
    clearplot(ds);
    updateplot(ds);
end

% Update the main fig axis limits to fit the new histograms
dfupdatexlim;
dfupdateylim;

if (v6 == true) % set default
    dfgetset('binDlgInfo', binDlgInfo);
end

