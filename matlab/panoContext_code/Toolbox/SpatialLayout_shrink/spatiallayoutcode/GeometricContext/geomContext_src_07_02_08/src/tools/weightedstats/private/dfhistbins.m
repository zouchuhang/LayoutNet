function [centers,edges] = dfhistbins(data,cens,freq,binInfo,F,x)
%DFHISTBINS Compute bin centers for a histogram
%   [CENTERS,EDGES] = DFHISTBINS(DATA,CENS,FREQ,BININFO,F,X) computes
%   histogram bin centers and edges for the rule specified in BININFO.  For
%   the Freedman-Diaconis rule, DFHISTBINS uses the empirical distribution
%   function F evaluated at the values X to compute the IQR.  When there is
%   censoring, DFHISTBINS cannot compute the Scott rule, and F-D is
%   substituted.

%   $Revision: 1.1.6.3 $  $Date: 2004/01/24 09:35:43 $
%   Copyright 2001-2004 The MathWorks, Inc.

xmin = min(data);
xmax = max(data);
xrange = xmax - xmin;
if isempty(freq)
    n = length(data);
else
    n = sum(freq);
end

rule = binInfo.rule;
% Can't compute the variance for the Scott rule when there is censoring,
% use F-D instead.
if (rule == 2) && ~isempty(cens) && any(cens)
    rule = 1; % Freedman-Diaconis
end

switch rule
case 1 % Freedman-Diaconis
    % Get "quartiles", which may not actually be the 25th and 75th points
    % if there is a great deal of censoring, and compute the IQR.
    iqr = diff(interp1q([F;1], [x;x(end)], [.25; .75]));
    
    % Guard against too small an IQR.  This may be because most
    % observations are censored, or because there are some extreme
    % outliers.
    if iqr < xrange ./ 10
        iqr = xrange ./ 10;
    end

    % Compute the bin width proposed by Freedman and Diaconis, and the
    % number of bins needed to span the data.  Use approximately that
    % many bins, placed at nice locations.
    [centers,edges] = binpicker(xmin, xmax, 'FD', n, iqr);

case 2 % Scott
    if isempty(freq)
        s = sqrt(var(data));
    else
        s = sqrt(var(data,freq));
    end

    % Compute the bin width proposed by Scott, and the number of bins
    % needed to span the data.  Use approximately that many bins,
    % placed at nice locations.
    [centers,edges] = binpicker(xmin, xmax, 'Scott', n, s);

case 3 % number of bins given
    % Do not create more than 1000 bins.
    [centers,edges] = binpicker(xmin, xmax, min(binInfo.nbins,1000));
    
case 4 % bins centered on integers
    xscale = max(abs([xmin xmax]));
    % If there'd be more than 1000 bins, center them on an appropriate
    % power of 10 instead.
    if xrange > 1000
        step = 10^ceil(log10(xrange/1000));
        xmin = step*round(xmin/step); % make the edges bin width multiples
        xmax = step*round(xmax/step);
        
    % If a bin width of 1 is effectively zero relative to the magnitude of
    % the endpoints, use a bigger power of 10.
    elseif xscale*eps > 1;
        step = 10^ceil(log10(xscale*eps));
        
    else
        step = 1;
    end
    centers = floor(xmin):step:ceil(xmax);
    edges = (floor(xmin)-.5*step):step:(ceil(xmax)+.5*step);
    
case 5 % bin width given
    % Do not create more than 1000 bins.
    binWidth = max(binInfo.width, xrange/1000);
    if (binInfo.placementRule == 1) % automatic placement: anchored at zero
        anchor = 0;
    else % anchored
        anchor = binInfo.anchor;
    end
    leftEdge = anchor + binWidth*floor((xmin-anchor) ./ binWidth);
    nbins = max(1,ceil((xmax-leftEdge) ./ binWidth));
    edges = leftEdge + (0:nbins) .* binWidth; % get exact multiples
    centers = edges(2:end) - 0.5 .* binWidth;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [centers,edges] = binpicker(xmin, xmax, nbins, nobs, extraArg)
%BINPICKER Generate pleasant bin locations for a histogram.
%   CENTERS = BINPICKER(XMIN,XMAX,NBINS) computes centers for histogram
%   bins spanning the range XMIN to XMAX, with extremes of the bins at
%   locations that are a multiple of 1, 2, 3, or 5 times a power of 10.
%
%   CENTERS = BINPICKER(XMIN,XMAX,'FD',N,IQR) uses the Freedman-Diaconis
%   rule for bin width to compute the number of bins.  N is the number of
%   data points, and IQR is the sample interquartile range of the data.
%
%   CENTERS = BINPICKER(XMIN,XMAX,'Scott',N,STD) uses Scott's rule for the
%   bin width to compute the number of bins.  N is the number of data
%   points, and STD is the sample standard deviation of the data.  Scott's
%   rule is appropriate for "normal-like" data.
%
%   CENTERS = BINPICKER(XMIN,XMAX,'Sturges',N) uses Sturges' rule for the
%   number of bins.  N is the number of data points.  Sturges' rule tends
%   to give fewer bins than either F-D or Scott.
%
%   For the Freedman-Diaconis, Scott's, or Sturges' rules, BINPICKER
%   automatically generates "nice" bin locations, where the bin width is 1,
%   2, 3, or 5 times a power of 10, and the bin edges fall on multiples of
%   the bin width.  Thus, the actual number of bins will often differ
%   somewhat from the number defined by the requested rule.
%
%   [CENTERS,EDGES] = BINPICKER(...) also returns the bin edges.

%   References:
%      [1] Freedman, D. and P. Diaconis (1981) "On the histogram as a
%          density estimator: L_2 theory", Zeitschrift fur
%          Wahrscheinlichkeitstheorie und verwandte Gebiete, 57:453–476.
%      [2] Scott, D.W. (1979) "On optimal and data-based histograms",
%          Biometrika, 66:605-610.
%      [3] Sturges, H.A. (1926) "The choice of a class interval",
%          J.Am.Stat.Assoc., 21:65-66.

if nargin < 3
    error('stats:binpicker:TooFewInputs', ...
          'Requires at least three inputs.');
elseif xmax < xmin
    error('stats:binpicker:MaxLessThanMin', ...
          'XMAX must be greater than or equal to XMIN.');
end

% Bin width rule specified
if ischar(nbins)
    ruleNames = ['fd     '; 'scott  '; 'sturges'];
    rule = strmatch(lower(nbins),ruleNames); % 1, 2, or 3
    if isempty(rule)
        error('stats:binpicker:UnknownRule', ...
              'RULE must be one of ''FD'', ''Scott'', or ''Sturges''.');
    elseif nobs < 1
        nbins = 1; % give 1 bin for zero-length data
        rule = 0;
    end

% Number of bins specified
else
    if nbins < 1 || round(nbins) ~= nbins
        error('stats:binpicker:NegativeNumBins', ...
            'NBINS must be a positive integer.');
    end
    rule = 0;
end

xscale = max(abs([xmin,xmax]));
xrange = xmax - xmin;

switch rule
case 1 % Freedman-Diaconis rule
    % Use the interquartile range to compute the bin width proposed by
    % Freedman and Diaconis, and the number of bins needed to span the
    % data.  Use approximately that many bins, placed at nice
    % locations.
    iqr = extraArg;
    rawBinWidth = 2*iqr ./ nobs.^(1/3);

case 2 % Scott's rule
    % Compute the bin width proposed by Scott, and the number of bins
    % needed to span the data.  Use approximately that many bins,
    % placed at nice locations.
    s = extraArg;
    rawBinWidth = 3.49*s ./ nobs.^(1/3);

case 3 % Sturges' rule for nbins
    nbins = 1 + log2(nobs);
    rawBinWidth = xrange ./ nbins;

otherwise % number of bins specified
    rawBinWidth = xrange ./ nbins;
end

% Make sure the bin width is not effectively zero.  Otherwise it will never
% amount to anything, which is what we knew all along.
rawBinWidth = max(rawBinWidth, eps*xscale);
% it may _still_ be zero, if data are all zeroes

% If the data are not constant, place the bins at "nice" locations
if xrange > max(sqrt(eps)*xscale, realmin)
    % Choose the bin width as a "nice" value.
    powOfTen = 10.^floor(log10(rawBinWidth)); % next lower power of 10
    relSize = rawBinWidth ./ powOfTen; % guaranteed in [1, 10)
    if  relSize < 1.5
        binWidth = 1*powOfTen;
    elseif relSize < 2.5
        binWidth = 2*powOfTen;
    elseif relSize < 4
        binWidth = 3*powOfTen;
    elseif relSize < 7.5
        binWidth = 5*powOfTen;
    else
        binWidth = 10*powOfTen;
    end

    % Automatic rule specified
    if rule > 0
        % Put the bin edges at multiples of the bin width, covering x.  The
        % actual number of bins used may not be exactly equal to the requested
        % rule. Always use at least two bins.
        leftEdge = binWidth*floor(xmin ./ binWidth);
        nbinsActual = max(2, ceil((xmax-leftEdge) ./ binWidth));

    % Number of bins specified
    else
        % Put the extreme bin edges at multiples of the bin width, covering x.
        % Then recompute the bin width to make the actual number of bins used
        % exactly equal to the requested number.
        leftEdge = binWidth*floor(xmin ./ binWidth);
        rightEdge = binWidth*ceil(xmax ./ binWidth);
        binWidth = (rightEdge - leftEdge) ./ nbins;
        nbinsActual = nbins;
    end

else % the data are nearly constant
    % For automatic rules, use a single bin.
    if rule > 0
        nbins = 1;
    end
    
    % There's no way to know what scale the caller has in mind, just create
    % something simple that covers the data.
    if xscale > realmin
        % Make the bins cover a unit width, or as small an integer width as
        % possible without the individual bin width being zero relative to
        % xscale.  Put the left edge on an integer or half integer below
        % xmin, with the data in the middle 50% of the bin.  Put the left
        % edge similarly above xmax.
        binRange = max(1, ceil(nbins*eps*xscale));
        leftEdge = floor(2*(xmin-binRange./4))/2;
        rightEdge = ceil(2*(xmax+binRange./4))/2;
    else
        leftEdge = -0.5;
        rightEdge = 0.5;
    end
    binWidth = (rightEdge - leftEdge) ./ nbins;
    nbinsActual = nbins;
end

edges = leftEdge + (0:nbinsActual) .* binWidth; % get exact multiples
centers = edges(2:end) - 0.5 .* binWidth;
