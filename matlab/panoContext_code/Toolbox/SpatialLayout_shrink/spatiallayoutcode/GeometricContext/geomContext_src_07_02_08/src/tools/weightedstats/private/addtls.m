function s = addtls(s)
%ADDTLS Add the t location-scale distribution.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.7 $  $Date: 2004/01/24 09:35:08 $

j = length(s) + 1;
s(j).name = 't location-scale';
s(j).code = 'tlocationscale';
s(j).pnames = {'mu' 'sigma' 'nu'};
s(j).pdescription = {'location' 'scale' 'shape'};
s(j).prequired = [false false false];
s(j).fitfunc = @tlsfit;
s(j).likefunc = @tlslike;
s(j).cdffunc = @tlscdf;
s(j).pdffunc = @tlspdf;
s(j).invfunc = @tlsinv;
s(j).statfunc = @tlsstat;
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = true;
s(j).paramvec = true;
s(j).support = [-Inf Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = false;
s(j).uselogpp = false;


% ==== t Location-Scale distribution functions ====

% these distribution functions do not yet handle arrays of parameters

function y = tlspdf(x, mu, sigma, nu)
%TLSPDF T location-scale probability density function (pdf).
sigma(sigma <= 0) = NaN;
nu(nu <= 0) = NaN;

y = tpdf((x - mu)./sigma,nu)./sigma;


function p = tlscdf(x ,mu, sigma, nu)
%TLSCDF T location-scale cumulative distribution function (cdf).
sigma(sigma <= 0) = NaN;
nu(nu <= 0) = NaN;

p = tcdf((x - mu)./sigma,nu);


function x = tlsinv(p, mu, sigma, nu)
%TLSINV Inverse of the t location-scale cumulative distribution function (cdf).
sigma(sigma <= 0) = NaN;
nu(nu <= 0) = NaN;

x = tinv(p,nu).*sigma + mu;


function r = tlsrnd(mu, sigma, nu, varargin)
%TLSRND Random arrays from the t location-scale distribution.
sigma(sigma <= 0) = NaN;
nu(nu <= 0) = NaN;

[err, sizeOut] = statsizechk(3,mu,sigma,nu,varargin{:});
if err > 0
    error('stats:tlsrnd:InconsistentSizes','Size information is inconsistent.');
end

r = mu + sigma.*trnd(nu,sizeOut);


function [m,v] = tlsstat(mu, sigma, nu)
%TLSSTAT Mean and variance for the t location-scale distribution.
sigma(sigma <= 0) = NaN;
nu(nu <= 0) = NaN;

if nu <= 1
    m = NaN;
else
    m = mu;
end
if nu <= 2
    v = Inf;
else
    v = sigma.^2 .* nu ./ (nu - 2);
end


function [nlogL,acov] = tlslike(params,data,cens,freq)
%TLSLIKE Negative log-likelihood for the t location-scale distribution.
if nargin < 4 || isempty(freq), freq = ones(size(data)); end
if nargin < 3 || isempty(cens), cens = zeros(size(data)); end

nlogL = tls_nloglf(params, data, cens, freq);
if nargout > 1
    acov = mlecov(params, data, 'nloglf',@tls_nloglf, 'cens',cens, 'freq',freq);
end


% ==== t location-scale fitting functions ====

function [phat,pci] = tlsfit(x,alpha,cens,freq,opts)

if nargin < 2 || isempty(alpha), alpha = .05; end
if nargin < 3 || isempty(cens), cens = zeros(size(x)); end
if nargin < 4 || isempty(freq), freq = ones(size(x)); end
if nargin < 5, opts = []; end

% Robust estimators for the mean and std dev of a normal, and method
% of moments on t-kurtosis for nu
xunc = x(cens == 0);
k = max(kurtosis(xunc), 4);
start = [median(xunc), 1.253.*mad(xunc), 2.*(2.*k-3)./(k-3)];

% The default options include turning fminsearch's display off.  This
% function gives its own warning/error messages, and the caller can turn
% display on to get the text output from fminsearch if desired.
options = statset(statset('tlsfit'), opts);
tolBnd = options.TolBnd;
options = optimset(options);

% Maximize the log-likelihood with respect to mu, sigma, and nu.
[phat,nll,err,output] = ...
    fminsearch(@tls_nloglf, start, options, x, cens, freq, tolBnd);
if (err == 0)
    % fminsearch may print its own output text; in any case give something
    % more statistical here, controllable via warning IDs.
    if output.funcCount >= options.MaxFunEvals
        wmsg = 'Maximum likelihood estimation did not converge.  Function evaluation limit exceeded.';
    else
        wmsg = 'Maximum likelihood estimation did not converge.  Iteration limit exceeded.';
    end
    if phat(3) > 100 % degrees of freedom became very large
       wmsg = sprintf('%s\n%s', wmsg, ...
                      'The normal distribution might provide a better fit.');
    end
    warning('stats:tlsfit:IterOrEvalLimit',wmsg);
elseif (err < 0)
    error('stats:tlsfit:NoSolution',...
          'Unable to reach a maximum likelihood solution.');
end

if nargout > 1
    acov = mlecov(phat, x, 'nloglf',@tls_nloglf, 'cens',cens, 'freq',freq);
    probs = [alpha/2; 1-alpha/2];
    se = sqrt(diag(acov))';

    % Compute the CI for mu using a normal approximation for muhat.
    pci(:,1) = norminv(probs, phat(1), se(1));

    % Compute the CI for sigma using a normal approximation for
    % log(sigmahat), and transform back to the original scale.
    % se(log(sigmahat)) is se(sigmahat) / sigmahat.
    logsigci = norminv(probs, log(phat(2)), se(2)./phat(2));
    pci(:,2) = exp(logsigci);

    % Compute the CI for nu using a normal distribution for nuhat.
    pci(:,3) = norminv(probs, phat(3), se(3));
end


function nll = tls_nloglf(parms, x, cens, freq, tolBnd)
%TLS_NLOGLF Objective function for t location-scale maximum likelihood.
mu = parms(1);
sigma = parms(2);
nu = parms(3);

% Restrict sigma and nu to the open interval (0, Inf).
if nargin > 4
    if sigma < tolBnd || nu < tolBnd
        nll = Inf;
        return
    end
end

t = (x - mu) ./ sigma;
w = nu + (t.^2);
logw = log(w);

L = -.5.*(nu+1).*logw + gammaln(.5.*(nu+1)) - gammaln(.5.*nu) + 0.5.*nu.*log(nu) - log(sigma) - .5.*log(pi);
ncen = sum(freq.*cens);
if ncen > 0
    cen = (cens == 1);
    if nu < 1e7  % Use the standard formula
        Scen = betainc(nu ./ w(cen), .5.*nu, 0.5) ./ 2;

        % Reflect for negative t.
        reflect = (t(cen) < 0);
        Scen(reflect) = 1 - Scen(reflect);

    else  % Use a normal approximation.
        Scen = log(0.5 * erfc(t(cen) ./ sqrt(2)));
    end
    L(cen) = log(Scen);
end
nll = -sum(freq .* L);

% Don't yet have dbetainc, so can't compute an analytic gradient with censoring.
%
% if nargout > 1
%     dL1 = (nu+1).*t./(w.*sigma);
%     dL2 = t.*dL1 - 1./sigma;
%     dL3 = .5.*(-logw - (nu+1)./w + psi(.5.*(nu+1)) - psi(.5.*nu) + log(nu) + 1);
%     if ncen > 0
% %         dL1(cen) = ;
% %         dL2(cen) = ;
% %         dL3(cen) = ;
%     end
%     ngrad = -[sum(freq .* dL1) sum(freq .* dL2) sum(freq .* dL3)];
% end
