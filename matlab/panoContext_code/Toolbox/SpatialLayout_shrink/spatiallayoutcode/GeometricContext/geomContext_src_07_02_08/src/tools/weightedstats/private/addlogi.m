function s = addlogi(s)
%ADDLOGI Add the logistic adistributions.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.9 $  $Date: 2004/01/24 09:35:06 $

j = length(s) + 1;
s(j).name = 'Logistic';
s(j).code = 'logistic';
s(j).pnames = {'mu' 'sigma'};
s(j).pdescription = {'location' 'scale'};
s(j).prequired = [false false];
s(j).fitfunc = @logifit;
s(j).likefunc = @logilike;
s(j).cdffunc = @logicdf;
s(j).pdffunc = @logipdf;
s(j).invfunc = @logiinv;
s(j).statfunc = @logistat;
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = true;
s(j).paramvec = true;
s(j).support = [-Inf Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = true;
s(j).uselogpp = false;

j = j + 1;
s(j).name = 'Log-Logistic';
s(j).code = 'loglogistic';
s(j).pnames = {'mu' 'sigma'};
s(j).pdescription = {'log location' 'log scale'};
s(j).prequired = [false false];
s(j).fitfunc = @loglfit;
s(j).likefunc = @logllike;
s(j).cdffunc = @loglcdf;
s(j).pdffunc = @loglpdf;
s(j).invfunc = @loglinv;
s(j).statfunc = @loglstat;
s(j).loginvfunc = @logiinv;
s(j).logcdffunc = @logicdf;
s(j).hasconfbounds = false;
s(j).censoring = true;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = true;
s(j).uselogpp = true;


% ==== Logistic distribution functions ====

% these distribution functions do not yet handle arrays of parameters

function y = logipdf(x, mu, sigma)
%LOGIPDF Logistic probability density function (pdf).
if (nargin<2), mu=0; end
if (nargin<3), sigma=1; end
sigma(sigma <= 0) = NaN;

z = (x - mu) ./ sigma;
k = (z>350); if any(k), z(k) = -z(k); end % prevent Inf/Inf
y = exp(z) ./ ((1 + exp(z)).^2 .* sigma);


function p = logicdf(x, mu, sigma)
%LOGICDF Logistic cumulative distribution function (cdf).
if (nargin<2), mu=0; end
if (nargin<3), sigma=1; end
sigma(sigma <= 0) = NaN;

p = 1 ./ (1 + exp(-(x - mu) ./ sigma));


function x = logiinv(p, mu, sigma)
%LOGIINV Inverse of the logistic cumulative distribution function (cdf).
if (nargin<2), mu=0; end
if (nargin<3), sigma=1; end
sigma(sigma <= 0) = NaN;

x = logit(p).*sigma + mu;


function r = logirnd(mu, sigma, varargin)
%LOGIRND Random arrays from the logistic distribution.
if (nargin<1), mu=0; end
if (nargin<2), sigma=1; end
sigma(sigma <= 0) = NaN;

[err, sizeOut] = statsizechk(2,mu,sigma,varargin{:});
if err > 0
    error('stats:logirnd:InconsistentSizes','Size information is inconsistent.');
end

p = rand(sizeOut);
r = log(p./(1-p)).*sigma + mu;


function [m,v] = logistat(mu, sigma)
%LOGISTAT Mean and variance for the logistic distribution.
if (nargin<1), mu=0; end
if (nargin<2), sigma=1; end
sigma(sigma <= 0) = NaN;

m = mu;
v = sigma.^2 .* pi.^2 ./ 3;


function [nlogL,acov] = logilike(params,data,cens,freq)
%LOGILIKE Negative log-likelihood for the logistic distribution.
if nargin < 4 || isempty(freq), freq = ones(size(data)); end
if nargin < 3 || isempty(cens), cens = zeros(size(data)); end

nlogL = logi_nloglf(params, data, cens, freq);
if nargout > 1
    acov = mlecov(params, data, 'nloglf',@logi_nloglf, 'cens',cens, 'freq',freq);
end


% ==== Logistic fitting functions ====

function [phat,pci] = logifit(x,alpha,cens,freq,opts)
%LOGIFIT Parameter estimates and confidence intervals for logistic data.

if nargin < 2 || isempty(alpha), alpha = .05; end
if nargin < 3 || isempty(cens), cens = zeros(size(x)); end
if nargin < 4 || isempty(freq), freq = ones(size(x)); end
if nargin < 5, opts = []; end

% Moment estimators as starting point
xunc = x(cens == 0);
start = [mean(xunc) std(xunc).*sqrt(3)./pi];

% The default options include turning statsfminbx's display off.  This
% function gives its own warning/error messages, and the caller can turn
% display on to get the text output from statsfminbx if desired.
options = statset(statset('logifit'), opts);
tolBnd = options.TolBnd;
options = optimset(options);
dfltOptions = struct('DerivativeCheck','off', 'HessMult',[], ...
    'HessPattern',ones(2,2), 'PrecondBandWidth',Inf, ...
    'TypicalX',ones(2,1), 'MaxPCGIter',1, 'TolPCG',0.1);

% Maximize the log-likelihood with respect to mu and sigma.
funfcn = {'fungrad' 'logifit' @logi_nloglf [] []};
[phat, nll, lagrange, err, output] = ...
         statsfminbx(funfcn, start, [-Inf; tolBnd], [Inf; Inf], ...
                     options, dfltOptions, 1, x, cens, freq);
if (err == 0)
    % statsfminbx may print its own output text; in any case give something
    % more statistical here, controllable via warning IDs.
    if output.funcCount >= options.MaxFunEvals
        wmsg = 'Maximum likelihood estimation did not converge.  Function evaluation limit exceeded.';
    else
        wmsg = 'Maximum likelihood estimation did not converge.  Iteration limit exceeded.';
    end
    warning('stats:logifit:IterOrEvalLimit',wmsg);
elseif (err < 0)
    error('stats:logifit:NoSolution',...
          'Unable to reach a maximum likelihood solution.');
end

if nargout > 1
    acov = mlecov(phat, x, 'nloglf',@logi_nloglf, 'cens',cens, 'freq',freq);
    probs = [alpha/2; 1-alpha/2];
    se = sqrt(diag(acov))';

    % Compute the CI for mu using a normal approximation for muhat.
    pci(:,1) = norminv(probs, phat(1), se(1));

    % Compute the CI for sigma using a normal approximation for
    % log(sigmahat), and transform back to the original scale.
    % se(log(sigmahat)) is se(sigmahat) / sigmahat.
    logsigci = norminv(probs, log(phat(2)), se(2)./phat(2));
    pci(:,2) = exp(logsigci);
end


function [nll,ngrad] = logi_nloglf(parms, x, cens, freq)
%LOGI_NLOGLF Objective function for logistic maximum likelihood.
mu = parms(1);
sigma = parms(2);
z = (x - mu) ./ sigma;
logitz = 1 ./ (1 + exp(-z));
clogitz = 1 ./ (1 + exp(z));
logclogitz = log(clogitz);
k = (z > 700); if any(k), logclogitz(k) = z(k); end % fix intermediate overflow

L = z + 2.*logclogitz - log(sigma);
ncen = sum(freq.*cens);
if ncen > 0
    cen = (cens == 1);
    L(cen) = logclogitz(cen);
end
nll = -sum(freq .* L);

if nargout > 1
    t = (2.*logitz - 1) ./ sigma;
    dL1 = t;
    dL2 = z.*t - 1./sigma;
    if ncen > 0
        t = logitz(cen) ./ sigma;
        dL1(cen) = t;
        dL2(cen) = z(cen) .* t;
    end
    ngrad = -[sum(freq .* dL1) sum(freq .* dL2)];
end



% ==== Log-Logistic distribution functions ====

% these distribution functions do not yet handle arrays of parameters

function y = loglpdf(x, mu, sigma)
%LOGLPDF Log-logistic probability density function (pdf).
if (nargin<2), mu=0; end
if (nargin<3), sigma=1; end
sigma(sigma <= 0) = NaN;

nonpos = (x <= 0);
x(nonpos) = realmin;
z = (log(x) - mu) ./ sigma;
c = ones(size(z));
k = (z>350); % prevent Inf/Inf
if any(k)
    z(k) = -z(k);
    c(k) = -1;
end
y = exp(z.*(1-c.*sigma) - mu) ./ ((1 + exp(z)).^2 .* sigma);
y(nonpos) = 0;
% the first and third of these would happen automatically for x==0, but
% generate LogOfZero warnings.  the second would be NaN.
y(x==0 & sigma<1) = 0;
y(x==0 & sigma==1) = 1;
y(x==0 & sigma>1) = Inf;


function p = loglcdf(x, mu, sigma)
%LOGLCDF Log-logistic cumulative distribution function (cdf).
if (nargin<2), mu=0; end
if (nargin<3), sigma=1; end
sigma(sigma <= 0) = NaN;

nonpos = (x <= 0);
x(nonpos) = realmin;
p = 1 ./ (1 + exp(-(log(x) - mu) ./ sigma));
% this would happen automatically for x==0, but generates LogOfZero warnings
p(nonpos) = 0;


function x = loglinv(p, mu, sigma)
%LOGLINV Inverse of the log-logistic cumulative distribution function (cdf).
if (nargin<2), mu=0; end
if (nargin<3), sigma=1; end
sigma(sigma <= 0) = NaN;

x = exp(logit(p).*sigma + mu);


function r = loglrnd(mu, sigma, varargin)
%LOGLRND Random arrays from the log-logistic distribution.
if (nargin<1), mu=0; end
if (nargin<2), sigma=1; end
sigma(sigma <= 0) = NaN;

[err, sizeOut] = statsizechk(2,mu,sigma,varargin{:});
if err > 0
    error('stats:loglrnd:InconsistentSizes','Size information is inconsistent.');
end

p = rand(sizeOut);
r = exp(log(p./(1-p)).*sigma + mu);


function [m,v] = loglstat(mu, sigma)
%LOGLSTAT Mean and variance for the log-logistic distribution.
if (nargin<1), mu=0; end
if (nargin<2), sigma=1; end
sigma(sigma <= 0) = NaN;

if sigma < 1
    m = exp(mu + gammaln(1+sigma) + gammaln(1-sigma));
else
    m = Inf;
end
if sigma < .5
    v = exp(2.*mu + gammaln(1+2.*sigma) + gammaln(1-2.*sigma)) - m.^2;
else
    v = Inf;
end


function [nlogL,acov] = logllike(params,data,cens,freq)
%LOGLLIKE Negative log-likelihood for the log-logistic distribution.
if nargin < 4 || isempty(freq), freq = ones(size(data)); end
if nargin < 3 || isempty(cens), cens = zeros(size(data)); end

nlogL = logl_nloglf(params, data, cens, freq);
if nargout > 1
    acov = mlecov(params, data, 'nloglf',@logl_nloglf, 'cens',cens, 'freq',freq);
end


% ==== Log-Logistic fitting functions ====

function [phat,pci] = loglfit(x,alpha,cens,freq,opts)
%LOGLFIT Parameter estimates and confidence intervals for log-logistic data.

if nargin < 2 || isempty(alpha), alpha = .05; end
if nargin < 3 || isempty(cens), cens = zeros(size(x)); end
if nargin < 4 || isempty(freq), freq = ones(size(x)); end
if nargin < 5, opts = []; end

if any(x <= 0)
    error('stats:loglfit:BadData','The data in X must be positive');
end

% Moment estimators as starting point
logxunc = log(x(cens == 0));
start = [mean(logxunc) std(logxunc).*sqrt(3)./pi];

% The default options include turning statsfminbx's display off.  This
% function gives its own warning/error messages, and the caller can turn
% display on to get the text output from statsfminbx if desired.
options = statset(statset('loglfit'), opts);
tolBnd = options.TolBnd;
options = optimset(options);
dfltOptions = struct('DerivativeCheck','off', 'HessMult',[], ...
    'HessPattern',ones(2,2), 'PrecondBandWidth',Inf, ...
    'TypicalX',ones(2,1), 'MaxPCGIter',1, 'TolPCG',0.1);

% Maximize the log-likelihood with respect to mu and sigma.
funfcn = {'fungrad' 'loglfit' @logl_nloglf [] []};
[phat, nll, lagrange, err, output] = ...
         statsfminbx(funfcn, start, [-Inf; tolBnd], [Inf; Inf], ...
                     options, dfltOptions, 1, x, cens, freq);
if (err == 0)
    % statsfminbx may print its own output text; in any case give something
    % more statistical here, controllable via warning IDs.
    if output.funcCount >= options.MaxFunEvals
        wmsg = 'Maximum likelihood estimation did not converge.  Function evaluation limit exceeded.';
    else
        wmsg = 'Maximum likelihood estimation did not converge.  Iteration limit exceeded.';
    end
    warning('stats:loglfit:IterOrEvalLimit',wmsg);
elseif (err < 0)
    error('stats:loglfit:NoSolution',...
          'Unable to reach a maximum likelihood solution.');
end

if nargout > 1
    acov = mlecov(phat, x, 'nloglf',@logi_nloglf, 'cens',cens, 'freq',freq);
    probs = [alpha/2; 1-alpha/2];
    se = sqrt(diag(acov))';

    % Compute the CI for mu using a normal approximation for muhat.
    pci(:,1) = norminv(probs, phat(1), se(1));

    % Compute the CI for sigma using a normal approximation for
    % log(sigmahat), and transform back to the original scale.
    % se(log(sigmahat)) is se(sigmahat) / sigmahat.
    logsigci = norminv(probs, log(phat(2)), se(2)./phat(2));
    pci(:,2) = exp(logsigci);
end


function [nll,ngrad] = logl_nloglf(parms, x, cens, freq)
%LOGL_NLOGLF Objective function for log-logistic maximum likelihood.
mu = parms(1);
sigma = parms(2);
logx = log(x);
z = (logx - mu) ./ sigma;
logitz = 1 ./ (1 + exp(-z));
clogitz = 1 ./ (1 + exp(z));
logclogitz = log(clogitz);
k = (z > 700); if any(k), logclogitz(k) = z(k); end % fix intermediate overflow

L = z + 2.*logclogitz - log(sigma) - logx;
ncen = sum(freq.*cens);
if ncen > 0
    cen = (cens == 1);
    L(cen) = logclogitz(cen);
end
nll = -sum(freq .* L);

if nargout > 1
    t = (2.*logitz - 1) ./ sigma;
    dL1 = t;
    dL2 = z.*t - 1./sigma;
    if ncen > 0
        t = logitz(cen) ./ sigma;
        dL1(cen) = t;
        dL2(cen) = z(cen) .* t;
    end
    ngrad = -[sum(freq .* dL1) sum(freq .* dL2)];
end


% ==== utility functions ====

function logitp = logit(p)
%LOGIT Logistic transformation, handling edge and out of range.
logitp = repmat(NaN,size(p));
logitp(p==0) = -Inf;
logitp(p==1) = Inf;
ok = (0<p & p<1);
logitp(ok) = log(p(ok)./(1-p(ok)));
