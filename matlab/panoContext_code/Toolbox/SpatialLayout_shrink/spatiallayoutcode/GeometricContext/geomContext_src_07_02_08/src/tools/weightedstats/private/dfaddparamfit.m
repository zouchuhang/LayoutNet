function hFit = dfaddparamfit(hFit, fitname, distname, dsname, fitframe, exclname, useestimated, fixedvals)
%DFADDPARAMFIT Add parametric fit in dfittool

%   $Revision: 1.1.6.10 $  $Date: 2004/01/24 09:35:10 $
%   Copyright 2003-2004 The MathWorks, Inc.

badfit = false;    % badfit=true means fit failed or not attempted
covok = true;      % covariance calculation can be done
if isempty(hFit)
    newfit = true;
    hFit = stats.dffit(fitname, fitframe);
else
    newfit = false;
end
listeners = hFit.listeners;
set(listeners, 'Enabled', 'off');

frameContents = fitframe.getContentPane;
componentsVector = get(frameContents,'Components');
fittingPanel = componentsVector(1);

% Get data set to fit
ds=find(getdsdb,'name',dsname);
hFit.distname = distname;
hFit.dataset  = dsname;
hFit.fittype = 'param'; 
hFit.dshandle = ds;

% Store some GUI values in fit
hFit.pfixedtext = fixedvals;
hFit.pestimated = useestimated;

% Extract data from this data set
alpha = 0.05;
hExcl = dfgetexclusionrule(exclname);
[x, cens, freq] = getincludeddata(ds,hExcl);

% Get information about the requested distribution
dist = dfgetdistributions(distname);
if length(dist)~=1 || isempty(x)
   if length(dist)~=1
      emsg = 'Bad distribution name.';
   else
      emsg = 'No data remaining after exclusion rule applied.';
   end
   wmsg = '';
   badfit = true;
end
if length(dist)==1
   hFit.enablebounds = dist.hasconfbounds;
end

% Perform the fit
lasterr('');
lastwarn('');
ws = warning('off');
if badfit
   p = [];
else
   try
      nparams = length(dist.pnames);
      if dist.censoring
         censargs = {cens freq};
      else
         if ~isempty(cens) && any(cens)
            error('stats:dfaddparamfit:NoCensoring',...
                  'Censoring not allowed with the %s distribution', distname);
         elseif ~isempty(freq) && any(freq~=1)
            x = expandInput(x,freq);
            freq = [];
         end
         censargs = {};
      end
   
      % How many output variables will this return?
      if dist.paramvec
         nparamvars = 1;
      else
         nparamvars = nparams;
      end
   
      fixedparams = cell(0,1);
      if any(~useestimated)
         for j=1:nparams
            if ~useestimated(j)
               txt = deblank(fixedvals{j});
               if isempty(txt)
                  error('stats:dfaddparamfit:BadParam',...
                        'Invalid value for parameter %s', dist.pnames{j});
               end
               num = str2double(txt);
               if ~isfinite(num)
                  error('stats:dfaddparamfit:BadParam',...
                        'Invalid value for parameter %s', dist.pnames{j});
               end
               fixedparams{length(fixedparams)+1} = num;
            end
         end
      end
   
      % Set up a cell array to receive outputs, then do the fit
      pcell = cell(nparamvars,1);
      [pcell{:}] = feval(dist.fitfunc, x, fixedparams{:}, alpha, censargs{:});

      % Extract results into a single vector
      if dist.paramvec
         p = pcell{1};
      else
         p = [pcell{:}];
      end
   catch
      p = [];
   end
end
warning(ws);

if ~badfit
   if ~isempty(lastwarn)
      wmsg = sprintf('Warning:  %s',lastwarn);
   else
      wmsg = '';
   end
   emsg = lasterr;
   newmsg = '';
   if any(~isfinite(p))
      newmsg = 'Fit produced infinite parameter estimates.';
   elseif numel(p)~=numel(dist.pnames) || ~isnumeric(p)
      newmsg = 'Fit function returned bad parameter values';
   end
   if ~isempty(newmsg)
      badfit = true;
      emsg = combinemsg(emsg,newmsg);
   end

   % Any type of failure so far makes the covariance calculation questionable
   if ~isempty(wmsg) || ~isempty(emsg)
      covok = false;
   end
end

% Try to get a likelihood value
if isempty(p)
   pcov = [];
   nlogl = NaN;
else
   try
      if ~isempty(dist.likefunc)
         if covok
            [nlogl,pcov] = feval(dist.likefunc, p, x, censargs{:});
         else
            nlogl = feval(dist.likefunc, p, x, censargs{:});
            pcov = [];
         end
      else
         pcov = [];
         nlogl = localnlogl(num2cell(p),dist.pdffunc,dist.cdffunc,x,cens,freq);
      end
      newmsg = '';
   catch
      newmsg = lasterr;
   end
   if isempty(newmsg) && (~isnumeric(nlogl) || ~isscalar(nlogl))
      newmsg = 'Result must be a numeric scalar';
   end
   if isnan(nlogl)
      nlogl = NaN;      % explicitly set to real nan to remove imaginary part
   end
   if ~isempty(newmsg);
      pcov = [];
      nlogl = NaN;
      wmsg = combinemsg(wmsg,...
                        sprintf('Error while evaluating likelihood:\n%s',...
                                newmsg));
   end
      
end

% Get the range over which to show the fit
dffig = dfgetset('dffig');
ax = findall(dffig,'Type','axes','Tag','main');
xlim = get(ax,'XLim');

% Create a fit object using the information we calculated
if badfit
   resultsText = emsg;
else
   try
      hFit = storefitresults(hFit, dist, p, pcov, nlogl, xlim, hExcl, exclname);
      resultsText = getresults(hFit);
   catch
      resultsText = lasterr;
      badfit = true;
   end
end

resultsText = combinemsg(wmsg,resultsText);

% Show results
hFit.resultstext = resultsText;
fittingPanel.setResults(resultsText)

if ~isempty(hFit)
   if ~newfit && ~(hFit.isgood == ~badfit)
   		com.mathworks.toolbox.stats.FitsManager.getFitsManager.fitIsGoodChanged(java(hFit), ~badfit);
   end
   hFit.isgood = ~badfit;
   if newfit
	  hFit.plot = 1;
      % Add to fit array
      connect(hFit,getfitdb,'up');
   end
end

if hFit.plot
   % Determine if bounds can be shown
   if ~dist.hasconfbounds
      hFit.showbounds = false;
   end
   
   % Update plotted curve
   updateplot(hFit);

   % Update plot limits
   dfswitchyard('dfupdatexlim');
   dfswitchyard('dfupdateylim');
end

set(listeners, 'Enabled', 'on');

if ~newfit
   com.mathworks.toolbox.stats.FitsManager.getFitsManager.fitChanged(...
       java(hFit),fitname,fitname);
end

% Display a more prominent warning outside the results text
if ~badfit && ~isempty(wmsg)
   warndlg(wmsg,'Distribution Fitting Warning','modal');
end

% ----------------------------------------------
function hFit = storefitresults(hFit, dist, p, pcov, nlogl, xlim, hExcl, exclname)
% Update its properties
hFit.distspec = dist;
hFit.params   = p;
hFit.pcov     = pcov;
hFit.pfixed   = false(size(p));
hFit.loglik = -nlogl;
hFit.support = dist.support;
hFit.exclusionrule = hExcl;
hFit.exclusionrulename = exclname;

hFit.xlim = xlim;
setftype(hFit,dfgetset('ftype'));


% ---------------------------------------------
function nlogl = localnlogl(p, pdf, cdf, x, cens, freq)
% Calculate negative log likelihood

% Handle defaults for option inputs
if isempty(cens)
   cens = false(size(x));
else
   cens = (cens == 1);
end
if isempty(freq)
   freq = ones(size(x));
end

% Compute for uncensored observations
nlogl = - sum(freq(~cens) .* log(feval(pdf, x(~cens), p{:})));

% Add component for censored observations
if any(cens)
   nlogl = nlogl - sum(freq(cens) .* log(1-feval(cdf, x(cens), p{:})));
end

% -----------------------------------
function msg = combinemsg(msg,newmsg)
%COMBINEMSG Combine multiple messages into a single message
if isempty(msg)
   msg = newmsg;
elseif ~isempty(newmsg)
   msg = sprintf('%s\n\n%s',msg,newmsg);
end

% -----------------------------------------
function expanded = expandInput(input,freq)
%EXPANDDATA Expand out an input vector using element frequencies.
if ~isequal(size(input),size(freq))
    error('stats:dfaddparamfit:InputSizeMismatch',...
          'Input argument sizes must match.');
end

% Remove points that have zero frequency
t = (freq == 0);
if any(t)
   input(t) = [];
   freq(t) = [];
end

% Expand the remainder
i = cumsum(freq);
j = zeros(1, i(end));
j(i(1:end-1)+1) = 1;
j(1) = 1;
expanded = input(cumsum(j));
