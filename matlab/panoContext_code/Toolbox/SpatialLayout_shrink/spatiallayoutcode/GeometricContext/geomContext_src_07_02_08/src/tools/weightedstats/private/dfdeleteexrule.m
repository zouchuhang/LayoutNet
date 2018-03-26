function dfdeleteexrule(names)
%DFDELETEEXRULE GUI helper to delete an exclusion rule

%   $Revision: 1.1.6.2 $  $Date: 2004/01/24 09:35:28 $
%   Copyright 2003-2004 The MathWorks, Inc.

fitdb = getfitdb;
fit = down(fitdb);
OKtoDelete = true;
fitsToDelete = {};

if ~isempty(fit)
    msg = '';
    for i=1:length(names)
        fitCnt = 0;
        fitnames = '';
        m='';
        fit = down(fitdb);
        while(~isempty(fit))
            if strcmp(names{i}, fit.exclusionrulename)
                fitsToDelete{1, end + 1} = fit.name;
                if fitCnt > 0
                    fitnames = [fitnames, ', '];
                end;
                fitCnt = fitCnt + 1;
                fitnames = [fitnames, fit.name];
            end
            fit = right(fit);
        end
        if fitCnt == 1
            m = sprintf('If you delete "%s", the following fit will also be deleted: %s\n', names{i}, fitnames);
        elseif fitCnt > 1
            m = sprintf('If you delete "%s", the following fits will also be deleted: %s\n', names{i}, fitnames);
        end
        msg = [msg, m]; 
    end
    if length(msg) > 0
        button = questdlg(msg, 'Deleting exclusion rules', 'OK', 'Cancel', 'OK');
        if ~strcmp(button, 'OK')
            OKtoDelete = false;
        end
    end 
end


if OKtoDelete
    import com.mathworks.toolbox.stats.*;
    if ~isempty(fitsToDelete)
        FitsManager.getFitsManager.deleteFits(fitsToDelete);
    end
    OutliersManager.getOutliersManager.deleteOutliers(names);
end





