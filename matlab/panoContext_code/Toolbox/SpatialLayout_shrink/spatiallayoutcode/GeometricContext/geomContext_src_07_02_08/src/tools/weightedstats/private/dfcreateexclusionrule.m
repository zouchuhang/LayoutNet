function err = dfcreateexclusionrule(name, dataset, yl, yh, yle, yhe)
%DFCREATEEXCLUSIONRULE Create dfittool exclusion rule

%   $Revision: 1.1.6.4 $  $Date: 2004/01/24 09:35:25 $
%   Copyright 2003-2004 The MathWorks, Inc.

err = '';
% check for duplicate name 
    outlierdb = getoutlierdb;
    outlier = down(outlierdb);

    while(~isempty(outlier))
        if strcmp(outlier.name, name)
            err = 'outliernamethesame';
            return;
        end
        outlier = right(outlier);
    end

% Make a new exclusion rule
stats.outlier(name, dataset, yl, yh, yle, yhe);



