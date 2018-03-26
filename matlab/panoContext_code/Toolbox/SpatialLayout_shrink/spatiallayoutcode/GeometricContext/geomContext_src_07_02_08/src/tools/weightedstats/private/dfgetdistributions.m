function [s,errid] = dfgetdistributions(distname,douser)
%DFGETDISTRIBUTIONS Get structure defining the distributions supported by dfittool

%   $Revision: 1.1.6.8 $  $Date: 2004/01/24 09:35:36 $
%   Copyright 2003-2004 The MathWorks, Inc.

errid = '';

% If a struct was passed in, store this for later use
if nargin>0 && isstruct(distname)
   dfgetset('alldistributions',distname);
   return
end

% Get old value if already created and stored
s = dfgetset('alldistributions');

% If not created yet, create it now
if isempty(s)
   % Get built-in distributions
   s = getbuiltins;

   if nargin<2 || douser
      % Get user-defined distributions (won't be done if we already
      % had a distribution list created before this function was called)
      [s,errid,errmsg] = dfgetuserdists(s);
      if ~isempty(errid)
         errordlg(errmsg,'DFITTOOL User-Defined Distributions','modal');
      end
   end

   % Sort by name
   lowernames = lower(strvcat(s.name));
   [ignore, ind] = sortrows(lowernames);
   s = s(ind);

   % Store it for next time
   dfgetset('alldistributions',s);
end

if nargin>0 && ~isempty(distname)
   % Return only the distribution(s) requested, not all of them
   allnames = {s.code};
   distnum = strmatch(lower(distname), allnames);
   s = s(distnum);
end


% ------------------------------------
function s = getbuiltins
%GETBUILTINS Get distributions functions provided off the shelf

ndists = 11;     % to be updated if distributions added or removed
s(ndists).name = '';

% Exponential distribution
j = 1;
s(j).name = 'Exponential';      % distribution name
s(j).code = 'exponential';      % distribution code name
s(j).pnames = {'mu'};           % parameter names
s(j).pdescription = {'scale'};  % parameter descriptions
s(j).prequired = false;         % is a value required for this parameter?
s(j).fitfunc = @expfit;         % fitting function
s(j).likefunc = @explike;       % likelihood (and covariance) function
s(j).cdffunc = @expcdf;         % cdf function
s(j).pdffunc = @exppdf;         % pdf function
s(j).invfunc = @expinv;         % inverse cdf function
s(j).statfunc = @expstat;       % function to compute mean and var
s(j).loginvfunc = [];           % inverse cdf function on log scale, if any
s(j).logcdffunc = [];           % cdf function on log scale, if any
s(j).hasconfbounds = true;      % supports conf bnds for cdf and inverse
s(j).censoring = true;          % supports censoring
s(j).paramvec = true;           % returns fitted parameters as a vector
s(j).support = [0 Inf];         % range of x with positive density
s(j).closedbound = [true false];% is x at this boundary point acceptable
s(j).iscontinuous = true;       % is continuous, not discrete
s(j).islocscale = true;         % is location/scale family, no shape param
s(j).uselogpp = false;          % use log scale for probability plot

% Extreme value
j = j+1;
s(j).name = 'Extreme value';
s(j).code = 'extreme value';
s(j).pnames = {'mu' 'sigma'};
s(j).pdescription = {'location' 'scale'};
s(j).prequired = [false false];
s(j).fitfunc = @evfit;
s(j).likefunc = @evlike;
s(j).cdffunc = @evcdf;
s(j).pdffunc = @evpdf;
s(j).invfunc = @evinv;
s(j).statfunc = @evstat;
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = true;
s(j).censoring = true;
s(j).paramvec = true;
s(j).support = [-Inf Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = true;
s(j).uselogpp = false;

% Gamma
j = j+1;
s(j).name = 'Gamma';
s(j).code = 'gamma';
s(j).pnames = {'a' 'b'};
s(j).pdescription = {'shape' 'scale'};
s(j).prequired = [false false];
s(j).fitfunc = @gamfit;
s(j).likefunc = @gamlike;
s(j).cdffunc = @gamcdf;
s(j).pdffunc = @gampdf;
s(j).invfunc = @gaminv;
s(j).statfunc = @gamstat;
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = true;
s(j).censoring = true;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = false;
s(j).uselogpp = false;

% Lognormal
j = j+1;
s(j).name = 'Lognormal';
s(j).code = 'lognormal';
s(j).pnames = {'mu' 'sigma'};
s(j).pdescription = {'log location' 'log scale'};
s(j).prequired = [false false];
s(j).fitfunc = @lognfit;
s(j).likefunc = @lognlike;
s(j).cdffunc = @logncdf;
s(j).pdffunc = @lognpdf;
s(j).invfunc = @logninv;
s(j).statfunc = @lognstat;
s(j).loginvfunc = @norminv;
s(j).logcdffunc = @normcdf;
s(j).hasconfbounds = true;
s(j).censoring = true;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = true;
s(j).uselogpp = true;

% Normal
j = j+1;
s(j).name = 'Normal';
s(j).code = 'normal';
s(j).pnames = {'mu' 'sigma'};
s(j).pdescription = {'location' 'scale'};
s(j).prequired = [false false];
s(j).fitfunc = @normfit;
s(j).likefunc = @normlike;
s(j).cdffunc = @normcdf;
s(j).pdffunc = @normpdf;
s(j).invfunc = @norminv;
s(j).statfunc = @normstat;
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = true;
s(j).censoring = true;
s(j).paramvec = false;
s(j).support = [-Inf Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = true;
s(j).uselogpp = false;

% Weibull
j = j+1;
s(j).name = 'Weibull';
s(j).code = 'weibull';
s(j).pnames = {'a' 'b'};
s(j).pdescription = {'scale' 'shape'};
s(j).prequired = [false false];
s(j).fitfunc = @wblfit;
s(j).likefunc = @wbllike;
s(j).cdffunc = @wblcdf;
s(j).pdffunc = @wblpdf;
s(j).invfunc = @wblinv;
s(j).statfunc = @wblstat;
s(j).loginvfunc = @evinv;
s(j).logcdffunc = @evcdf;
s(j).hasconfbounds = true;
s(j).censoring = true;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = true;
s(j).uselogpp = true;

% Rayleigh
j = j+1;
s(j).name = 'Rayleigh';
s(j).code = 'rayleigh';
s(j).pnames = {'b'};
s(j).pdescription = {'scale'};
s(j).prequired = false;
s(j).fitfunc = @raylfit;
s(j).likefunc = [];
s(j).cdffunc = @raylcdf;
s(j).pdffunc = @raylpdf;
s(j).invfunc = @raylinv;
s(j).statfunc = @raylstat;
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = false;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = true;
s(j).uselogpp = false;

% Poisson
j = j+1;
s(j).name = 'Poisson';
s(j).code = 'poisson';
s(j).pnames = {'lambda'};
s(j).pdescription = {'mean'};
s(j).prequired = false;
s(j).fitfunc = @poissfit;
s(j).likefunc = [];
s(j).cdffunc = @poisscdf;
s(j).pdffunc = @poisspdf;
s(j).invfunc = @poissinv;
s(j).statfunc = @poisstat;
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = false;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [true false];
s(j).iscontinuous = false;
s(j).islocscale = false;
s(j).uselogpp = false;

% Negative binomial
j = j+1;
s(j).name = 'Negative Binomial';
s(j).code = 'negative binomial';
s(j).pnames = {'r' 'p'};
s(j).pdescription = {'' ''};
s(j).prequired = [false false];
s(j).fitfunc = @nbinfit;
s(j).likefunc = @nbinlike;
s(j).cdffunc = @nbincdf;
s(j).pdffunc = @nbinpdf;
s(j).invfunc = @nbininv;
s(j).statfunc = @nbinstat;
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = false;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [true false];
s(j).iscontinuous = false;
s(j).islocscale = false;
s(j).uselogpp = false;

% Beta
j = j+1;
s(j).name = 'Beta';
s(j).code = 'beta';
s(j).pnames = {'a' 'b'};
s(j).pdescription = {'' ''};
s(j).prequired = [false false];
s(j).fitfunc = @betafit;
s(j).likefunc = @betalike;
s(j).cdffunc = @betacdf;
s(j).pdffunc = @betapdf;
s(j).invfunc = @betainv;
s(j).statfunc = @betastat;
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = false;
s(j).paramvec = true;
s(j).support = [0 1];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = false;
s(j).uselogpp = false;

% Binomial
j = j+1;
s(j).name = 'Binomial';
s(j).code = 'binomial';
s(j).pnames = {'N' 'p'};
s(j).pdescription = {'trials' 'probability'};
s(j).prequired = [true false];
s(j).fitfunc = @localbinofit;
s(j).likefunc = [];
s(j).cdffunc = @binocdf;
s(j).pdffunc = @binopdf;
s(j).invfunc = @binoinv;
s(j).statfunc = @binostat;
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = false;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [true false];
s(j).iscontinuous = false;
s(j).islocscale = false;
s(j).uselogpp = false;

s = addbisa(s);
s = addinvg(s);
s = addlogi(s);
s = addnaka(s);
s = addtls(s);
s = addrice(s);

% ------------ binomial function is a special case
function [phat,pci] = localbinofit(x,N,alpha)
%LOCALBINOFIT Version of binofit that operates on vectors

nx = length(x);
sumx = sum(x);
sumN = nx * N;
if nargout==2
   [phat,pci] = binofit(sumx,sumN,alpha);
else
   phat = binofit(sumx,sumN,alpha);
end

phat = [N phat];
if nargout==2
   pci = [NaN NaN; pci]';
end
