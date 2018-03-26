function [errmsg,x,values] = dfevaluate(fitNames,x,fun,wantBounds,confLevel,plotFun,dum)
%DFEVALUATE Evaluate fits for DFITTOOL

%   $Revision: 1.1.6.3 $  $Date: 2004/01/24 09:35:32 $
%   Copyright 1993-2004 The MathWorks, Inc.


% If the function is th empty string, clear the plot (if there is one) and
% clear any saved data.
if isempty(fitNames)
    dfevaluateplot(false,false); % closes the plot window
    dfgetset('evaluateResults', []);
    return
end
    
nfits = length(fitNames);
try
    x = sprintf('[ %s ]', x); % allow an unbracketed list of numbers to work
    x = eval(x);
    confLevel = eval(confLevel) ./ 100;
catch
    x = [];
    values = zeros(0, nfits*(1+2*wantBounds));
    errmsg = sprintf('Invalid MATLAB expression: %s',lasterr);
    h = [];
    return
end
errmsg = '';

x = x(:);
n = length(x);

% % This is enforced by the evaluate panel.
% switch fun
% case {'pdf' 'hazrate' 'condmean'}
%     % No bounds allowed for pdf, hazrate, or conditional mean.
%     wantBounds = false;
% otherwise % {'cdf' 'icdf' 'survivor' 'cumhazard' 'probplot'}
%     % bounds are allowed
% end   

% Output table will have first column for data, then for each fit, one
% column for function, two columns for bounds (if requested).
values = repmat(NaN, n, nfits*(1+2*wantBounds));

fitdb = getfitdb;
for i = 1:nfits
    fit = find(fitdb, 'name', fitNames{i});
    % Cannot compute bounds for kernel smooths and certain parametric fits.
    getBounds = wantBounds && ...
        (~isequal(fit.fittype, 'smooth') && fit.distspec.hasconfbounds);

    % Evaluate the requested function for this fit.
    if getBounds
        [y,ylo,yup] = eval(fit,x,fun,confLevel);
        values(:,3*i-2) = y;
        values(:,3*i-1) = ylo;
        values(:,3*i) = yup;
    else
        y = eval(fit,x,fun);
        if wantBounds
            values(:,3*i-2) = y;
        else
            values(:,i) = y;
        end
    end
end

% Save the results for SaveToWorkSpace.  The plot function can be called
% directly from java, so it uses those saved results as well.
dfgetset('evaluateResults', [x,values]);

% Save information about the fits that we've evaluated for the plot function.
dfgetset('evaluateFun', fun);
dfgetset('evaluateInfo', struct('fitNames',{fitNames}, 'wantBounds',wantBounds));

dfevaluateplot(plotFun,dum);
