function [imsource, D, C, F, yexp, censexp, freqexp] = dfviewdatapreview(dataset, width, height)
% DFVIEWDATA Helper function for the dfittool viewdata panel
%
%    [X, Y, W] = DFVIEWDATA(DATASET)
%    returns the x, y and w values for the given dataset
%	 (in a manner that the Java GUI can use)

%   Copyright 2003-2004 The MathWorks, Inc. 
%   $Revision: 1.1.6.5 $  $Date: 2004/01/24 09:36:06 $

ds = handle(dataset);

[err, imsource] = dfpreview('', '', '', width, height, ds);
                            
[D, C, F] = dfviewdata(dataset);
C = double(C); %C might be logical - convert to double for GUI compatibility 

yexp = ds.yexp;
censexp = ds.censexp;
freqexp = ds.freqexp;

