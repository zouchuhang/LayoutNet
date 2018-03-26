function dfevaluateplot(plotFun,dum)
%DFEVALUATEPLOT Plot data and evaluated fits for DFITTOOL

%   $Revision: 1.1.6.3 $  $Date: 2004/01/24 09:35:33 $
%   Copyright 1993-2004 The MathWorks, Inc.

plotfig = dfgetset('evaluateFigure');

% If no plotting selected, delete the existing figure if there is one
if ~plotFun % && ~plotData
    if ~isempty(plotfig) && ishandle(plotfig)
%         h = findobj(allchild(plotfig),'flat','serializable','on');
%         delete(h);
%         plotaxes = axes('Visible','on', 'XTick',[], 'YTick',[], 'Parent',plotfig);
%         text(.5, .5, xlate('Press "Apply" to create a new plot'), ...
%              'Parent',plotaxes, 'HorizontalAlignment','center');
        delete(plotfig);
        dfgetset('evaluateFigure',plotfig);
    end
    return;
end

% Get the current evaluated results
evaluateResults = dfgetset('evaluateResults');
x = evaluateResults(:,1);
values = evaluateResults(:,2:end);

% Get the current information about the evaluated results
evaluateInfo = dfgetset('evaluateInfo');
fitNames = evaluateInfo.fitNames;
wantBounds = evaluateInfo.wantBounds;
fun = dfgetset('evaluateFun');

nfits = length(fitNames);

% Create plotting figure if it does not yet exist
if isempty(plotfig) || ~ishandle(plotfig)
    plotfig = figure('Visible','on', ...
                     'IntegerHandle','off',...
                     'HandleVisibility','callback',...
                     'name','Distribution Fitting Evaluate',...
                     'numbertitle','off',...
                     'PaperPositionMode','auto',...
                     'doublebuffer','on',...
                     'CloseRequestFcn',@closefig);
    dfgetset('evaluateFigure',plotfig);
end

% New or old, prepare figure by removing old contents
h = findobj(allchild(plotfig),'flat','serializable','on');
delete(h);

plotaxes = axes('Parent',plotfig);

% Will save fit line handles for legend, but not conf bound handles.
lineHndls = repmat(NaN,nfits,1);

% % Need to save fit/dataset names and line handles for legend
% nhndls = nfits * (plotFun + plotData); % might multiple count some datasets
% lineHndls = repmat(NaN,nhndls,1);
% lgndNames = cell(nhndls,1);
% lgndOrder = zeros(nhndls,1);
% nfitHndls = nfits * plotFun; % how many fits will be plotted
% lineCnt = nfitHndls; % keep track of how many (unique) legend items

% If there's only one point to plot, make it more visible.
if isscalar(x)
    marker = '.';
else
    marker = 'none';
end

fitdb = getfitdb;
for i = 1:nfits
    fit = find(fitdb, 'name', fitNames{i});
    getBounds = wantBounds && ...
        (~isequal(fit.fittype, 'smooth') && fit.distspec.hasconfbounds);

    if wantBounds
        y = values(:,3*i-2);
        if getBounds
            ylo = values(:,3*i-1);
            yup = values(:,3*i);
        end
    else
        y = values(:,i);
    end

    % Plot the function (and bounds) for this fit.
    if plotFun
%         lgndNames{i} = fit.name;
        color = fit.ColorMarkerLine{1};
        lineHndls(i) = line(x,y, 'LineStyle','-', ...
            'Marker',marker, 'Color',color, 'Parent',plotaxes);
        if getBounds
            line(x,ylo, 'LineStyle','--', ...
                'Marker',marker, 'Color',color, 'Parent',plotaxes);
            line(x,yup, 'LineStyle','--', ...
                'Marker',marker, 'Color',color, 'Parent',plotaxes);
        end

%         % If data are not being plotted, the fits will appear in the legend
%         % in their natural order from fitNames.  Otherwise, lgndOrder will
%         % be modified below so that they will appear under their dataset.
%         lgndOrder(i) = i;
    end

%     % Plot the empirical function (and bounds) for the data.
%     if plotData
%         % First figure out if this dataset has already been plotted
%         ds = fit.dshandle;
%         dataIdx = strmatch(ds.name, lgndNames((nfitHndls+1):lineCnt),'exact');
%         if isempty(dataIdx)
%             lineCnt = lineCnt + 1; % one more line (dataset) has been plotted
%             dataIdx = lineCnt - nfitHndls;
%             lgndNames{lineCnt} = ds.name;
%             color = ds.ColorMarkerLine{1};
%             if wantBounds && ~isequal(fun,'icdf')
%                 % *** this uses the current conf level saved in the dataset, not
%                 % *** the one requested from the evaluate panel
%                 [xdata,ydata,ydatabounds] = getplotdata(ds,fun,false);
%                 lineHndls(lineCnt) = line(xdata,ydata, 'LineStyle','-', ...
%                     'Marker','none', 'Color',color, 'Parent',plotaxes);
%                 line(xdata,ydatabounds(:,1), 'LineStyle','--', ...
%                     'Marker','none', 'Color',color, 'Parent',plotaxes);
%                 line(xdata,ydatabounds(:,2), 'LineStyle','--', ...
%                     'Marker','none', 'Color',color, 'Parent',plotaxes);
%             else
%                 [xdata,ydata] = getplotdata(ds,fun);
%                 lineHndls(lineCnt) = line(xdata,ydata, 'LineStyle','-', ...
%                     'Marker','none', 'Color',color, 'Parent',plotaxes);
%             end
%
%             % Datasets will appear in the legend in their natural order
%             % from fitNames, and the corresponding fits will appear
%             % directly below.  If more than one fit use the same dataset,
%             % that dataset will appear only once.
%             lgndOrder(lineCnt) = 10000*dataIdx;
%         end
%
%         % Put the current fit below its dataset in the legend.
%         if plotData
%             lgndOrder(i) = 10000*dataIdx + i;
%         end
%     end
end

if plotFun % || plotData
%     set(plotaxes, 'XLim', [min(x), max(x)]);
%     [dum,ord] = sort(lgndOrder(1:lineCnt));
%     legend(plotaxes,lineHndls(ord),lgndNames(ord),0);
    legend(plotaxes,lineHndls,fitNames,0);
    figure(plotfig);
end

% ------------------------------
function closefig(varargin)
%CLOSEFIG Close this figure, but also update check box in evaluate panel

% Delete the figure containing the evaluate plot
h = gcbf;
if ~isempty(h) && ishandle(h)
   delete(h);
end

% Update the checkbox
evp = com.mathworks.toolbox.stats.Evaluate.getEvaluatePanel;
if ~isempty(evp) && ~isequal(evp,0)
   evp.setPlotCB(false);
end
