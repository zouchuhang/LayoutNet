function [v0, v1, v2, v3, v4, v5, v6, imsource] = dbinwidthpreview(dataset, width, height)
% DBINWIDTHPREVIEWDATA Helper function for the dfittool set bin width panel
%

%   Copyright 2003-2004 The MathWorks, Inc. 
%   $Revision: 1.1.6.2 $  $Date: 2004/01/24 09:35:17 $

ds = handle(dataset);

[err, imsource] = dfpreview('', '', '', width, height, ds);
binDlgInfo = ds.binDlgInfo;
v0 = binDlgInfo.rule;
v1 = binDlgInfo.nbinsExpr;
v2 = binDlgInfo.widthExpr;
v3 = binDlgInfo.placementRule;
v4 = binDlgInfo.anchorExpr;
v5 = binDlgInfo.applyToAll;
v6 = binDlgInfo.setDefault;


