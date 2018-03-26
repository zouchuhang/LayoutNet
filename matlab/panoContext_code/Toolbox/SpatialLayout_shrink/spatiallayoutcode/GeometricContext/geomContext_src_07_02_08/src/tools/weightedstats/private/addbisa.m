function s = addbisa(s)
%ADDBISA Add the Birnbaum-Saunders distribution.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.8 $  $Date: 2004/01/24 09:35:04 $

j = length(s) + 1;
s(j).name = 'Birnbaum-Saunders';
s(j).code = 'birnbaumsaunders';
s(j).pnames = {'beta' 'gamma'};
s(j).pdescription = {'scale' 'shape'};
s(j).prequired = [false false];
s(j).fitfunc = @bisafit;
s(j).likefunc = @bisalike;
s(j).cdffunc = @bisacdf;
s(j).pdffunc = @bisapdf;
s(j).invfunc = @bisainv;
s(j).statfunc = @bisastat;
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = true;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = false;
s(j).uselogpp = false;


% ==== Birnbaum-Saunders distribution functions ====

% these distribution functions do not yet handle arrays of parameters

function y = bisapdf(x, beta, gamma)
%BISAPDF Birnbaum-Saunders probability density function (pdf).
beta(beta <= 0) = NaN;
gamma(gamma <= 0) = NaN;

nonpos = (x <= 0);
x(nonpos) = realmin;
z = (sqrt(x./beta) - sqrt(beta./x)) ./ gamma;
w = (sqrt(x./beta) + sqrt(beta./x)) ./ gamma;
ynorm = exp(-0.5 .* z.^2) ./ sqrt(2.*pi);
y = ynorm .* w ./ (2.*x);
% this would happen automatically for x==0, but generates DivideByZero warnings
y(nonpos) = 0;


function p = bisacdf(x, beta, gamma)
%BISACDF Birnbaum-Saunders cumulative distribution function (cdf).
beta(beta <= 0) = NaN;
gamma(gamma <= 0) = NaN;

nonpos = (x <= 0);
x(nonpos) = realmin;
z = (sqrt(x./beta) - sqrt(beta./x)) ./ gamma;
p = 0.5 * erfc(-z ./ sqrt(2));
% this would happen automatically for x==0, but generates DivideByZero warnings
p(nonpos) = 0;


function x = bisainv(p, beta, gamma)
%BISAINV Inverse of the Birnbaum-Saunders cumulative distribution function (cdf).
beta(beta <= 0) = NaN;
gamma(gamma <= 0) = NaN;

p(p < 0 | 1 < p) = NaN;
gamz = -sqrt(2).*erfcinv(2*p) .* gamma;
x = 0.25 .* beta .* (gamz + sqrt(4+gamz.^2)).^2;
x(p == 0 & ~isnan(beta) & ~isnan(gamma)) = 0;
% x(p == 1 & ~isnan(beta) & ~isnan(gamma)) = Inf; % get this automatically

function r = bisarnd(beta, gamma, varargin)
%BISARND Random arrays from the Birnbaum-Saunders distribution.
beta(beta <= 0) = NaN;
gamma(gamma <= 0) = NaN;

[err, sizeOut] = statsizechk(2,beta,gamma,varargin{:});
if err > 0
    error('stats:bisarnd:InconsistentSizes','Size information is inconsistent.');
end

plusminus = 2.*(rand(sizeOut)>.5) - 1; % plus or minus one, w.p. 1/2
gamz = gamma.*randn(sizeOut);
r = 0.5.*beta .* (2 + gamz.^2 + plusminus.*gamz.*sqrt(4 + gamz.^2));
% gamz = randn(sizeOut) .* gamma;
% r = 0.25 .* beta .* (gamz + sqrt(4+gamz.^2)).^2;


function [m,v] = bisastat(beta, gamma)
%BISASTAT Mean and variance for the Birnbaum-Saunders distribution.
beta(beta <= 0) = NaN;
gamma(gamma <= 0) = NaN;

m = beta .* (0.5 .* gamma.^2 + 1);
v = (beta.*gamma).^2 .* (1.25 .* gamma.^2 + 1);


function [nlogL,acov] = bisalike(params,data,cens,freq)
%BISALIKE Negative log-likelihood for the Birnbaum-Saunders distribution.
if nargin < 4 || isempty(freq), freq = ones(size(data)); end
if nargin < 3 || isempty(cens), cens = zeros(size(data)); end

nlogL = bisa_nloglf(params, data, cens, freq);
if nargout > 1
    acov = mlecov(params, data, 'nloglf',@bisa_nloglf, 'cens',cens, 'freq',freq);
end


% ==== Birnbaum-Saunders fitting functions ====

function [phat,pci] = bisafit(x,alpha,cens,freq,opts)
%BISAFIT Parameter estimates and confidence intervals for Birnbaum-Saunders data.

if nargin < 2 || isempty(alpha), alpha = .05; end
if nargin < 3 || isempty(cens), cens = zeros(size(x)); end
if nargin < 4 || isempty(freq), freq = ones(size(x)); end
if nargin < 5, opts = []; end

if any(x <= 0)
    error('stats:bisafit:BadData','The data in X must be positive');
end

% Starting points as suggested by Birnbaum and Saunders
xunc = x(cens==0); xbarunc = mean(xunc); xinvbarunc = mean(1./xunc);
start = [sqrt(xbarunc./xinvbarunc) 2.*sqrt(sqrt(xbarunc.*xinvbarunc) - 1)];

% The default options include turning statsfminbx's display off.  This
% function gives its own warning/error messages, and the caller can turn
% display on to get the text output from statsfminbx if desired.
options = statset(statset('bisafit'), opts);
tolBnd = options.TolBnd;
options = optimset(options);
dfltOptions = struct('DerivativeCheck','off', 'HessMult',[], ...
    'HessPattern',ones(2,2), 'PrecondBandWidth',Inf, ...
    'TypicalX',ones(2,1), 'MaxPCGIter',1, 'TolPCG',0.1);

% Maximize the log-likelihood with respect to mu and sigma.
funfcn = {'fungrad' 'bisafit' @bisa_nloglf [] []};
[phat, nll, lagrange, err, output] = ...
         statsfminbx(funfcn, start, [tolBnd; tolBnd], [Inf; Inf], ...
                     options, dfltOptions, 1, x, cens, freq);
if (err == 0)
    % statsfminbx may print its own output text; in any case give something
    % more statistical here, controllable via warning IDs.
    if output.funcCount >= options.MaxFunEvals
        wmsg = 'Maximum likelihood estimation did not converge.  Function evaluation limit exceeded.';
    else
        wmsg = 'Maximum likelihood estimation did not converge.  Iteration limit exceeded.';
    end
    warning('stats:bisafit:IterOrEvalLimit',wmsg);
elseif (err < 0)
    error('stats:bisafit:NoSolution',...
          'Unable to reach a maximum likelihood solution.');
end

% Compute CIs using a normal approximation for phat.
if nargout > 1
    acov = mlecov(phat, x, 'nloglf',@bisa_nloglf, 'cens',cens, 'freq',freq);
    probs = [alpha/2; 1-alpha/2];
    se = sqrt(diag(acov))';
    pci = norminv([probs probs], [phat; phat], [se; se]);
end


function [nll,ngrad] = bisa_nloglf(params, data, cens, freq)
%BISA_NLOGLF Objective function for Birnbaum-Saunders maximum likelihood.

beta = params(1);
gamma = params(2);
z = (sqrt(data./beta) - sqrt(beta./data)) ./ gamma;
w = (sqrt(data./beta) + sqrt(beta./data)) ./ gamma;

logphi = -0.5 .* (z.^2 + log(2.*pi));
L = logphi + log(w) - log(2.*data);
ncen = sum(freq.*cens);
if ncen > 0
    cen = (cens == 1);
    zcen = z(cen);
    Scen = 0.5 * erfc(zcen ./ sqrt(2));
    L(cen) = log(Scen);
end
nll = -sum(freq .* L);

if nargout > 1
    dL1 = (w.^2 - 1) .* 0.5.*z./(w.*beta);
    dL2 = (z.^2 - 1) ./ gamma;
    if ncen > 0
        phicen = exp(logphi(cen));
        wcen = w(cen);
        d1Scen = phicen .* 0.5.*wcen./beta;
        d2Scen = phicen .* zcen./gamma;
        dL1(cen) = d1Scen ./ Scen;
        dL2(cen) = d2Scen ./ Scen;
    end
    ngrad = -[sum(freq .* dL1) sum(freq .* dL2)];
end
