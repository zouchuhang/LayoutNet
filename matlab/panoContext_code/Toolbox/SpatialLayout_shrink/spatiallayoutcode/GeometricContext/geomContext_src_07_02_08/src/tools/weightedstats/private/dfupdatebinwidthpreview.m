function [err, imsource] = dfupdatebinwidthpreview(dataset, width, height, v0, v1, v2, v3, v4)
% DFUPDATEBINWIDTHPREVIEW Helper function for dfittool

%   Copyright 2003-2004 The MathWorks, Inc. 
%   $Revision: 1.1.6.2 $  $Date: 2004/01/24 09:36:00 $

ds = handle(dataset);
rules.rule = v0;
rules.nbins = str2num(v1);
rules.width = str2num(v2);
rules.placementRule = v3;
rules.anchor = str2num(v4);

[err, imsource] = dfpreview('', '', '', width, height, ds, rules);

