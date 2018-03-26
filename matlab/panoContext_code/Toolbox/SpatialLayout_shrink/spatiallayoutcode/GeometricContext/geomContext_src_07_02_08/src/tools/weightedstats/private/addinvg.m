function s = addinvg(s)
%ADDINVG Add the inverse Gaussian distribution.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.10 $  $Date: 2004/04/15 01:01:53 $

j = length(s) + 1;
s(j).name = 'Inverse Gaussian';
s(j).code = 'inversegaussian';
s(j).pnames = {'mu' 'lambda'};
s(j).pdescription = {'scale' 'shape'};
s(j).prequired = [false false];
s(j).fitfunc = @invgfit;
s(j).likefunc = @invglike;
s(j).cdffunc = @invgcdf;
s(j).pdffunc = @invgpdf;
s(j).invfunc = @invginv;
s(j).statfunc = @invgstat;
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


% ==== inverse Gaussian distribution functions ====

% these distribution functions do not yet handle arrays of parameters

function y = invgpdf(x, mu, lambda)
%INVGPDF inverse Gaussian probability density function (pdf).
mu(mu <= 0) = NaN;
lambda(lambda <= 0) = NaN;

nonpos = (x <= 0);
x(nonpos)= realmin;
y = sqrt(lambda./(2.*pi.*x.^3)) .* exp(-0.5.*lambda.*(x./mu - 2 + mu./x)./mu);
% this would happen automatically for x==0, but generates DivideByZero warnings
y(nonpos) = 0;


function p = invgcdf(x, mu, lambda)
%INVGCDF inverse Gaussian cumulative distribution function (cdf).
mu(mu <= 0) = NaN;
lambda(lambda <= 0) = NaN;

nonpos = (x <= 0);
x(nonpos)= realmin;
z1 = (x./mu - 1).*sqrt(lambda./x);
z2 = -(x./mu + 1).*sqrt(lambda./x);
p = 0.5.*erfc(-z1./sqrt(2)) + exp(2.*lambda./mu) .* 0.5.*erfc(-z2./sqrt(2));
% this would happen automatically for x==0, but generates DivideByZero warnings
p(nonpos) = 0;


function x = invginv(p, mu, lambda)
%INVGINV Inverse of the inverse Gaussian cumulative distribution function (cdf).
mu(mu <= 0) = NaN;
lambda(lambda <= 0) = NaN;

k = (0 < p & p < 1);
allOK = all(k(:));
if isa(p,'single')
   h = ones(size(p),'single');
   reltol = eps('single').^(3/4);
   mynan = single(NaN);
else
   h = ones(size(p));
   reltol = eps.^(3/4);
   mynan = NaN;
end

% Fill in NaNs for out of range cases, fill in edges cases when P is 0 or 1.
if ~allOK
    x = repmat(mynan, size(k));
    x(p == 0 & ~isnan(mu) & ~isnan(lambda)) = 0;
    x(p == 1 & ~isnan(mu) & ~isnan(lambda)) = Inf;
    
    % Remove the bad/edge cases, leaving the easy cases.  If there's
    % nothing remaining, return.
    if any(k(:))
        if numel(p) > 1
           p = p(k);
           h = h(k);
        end
    else
        return;
    end
end
    
% Newton's Method to find a root of invgcdf(x,mu,lambda) = p
%
% Choose a starting guess for q.  Use quantiles from a lognormal
% distribution with the same mean (==1) and variance (==lambda0) as
% IG(1,lambda0).
lambda0 = lambda/mu;
sigsqLN = log(1./lambda0 + 1);
muLN = -0.5 .* sigsqLN;
q = exp(muLN - sqrt(2.*sigsqLN).*erfcinv(2*p));

% Break out of the iteration loop when the relative size of the last step
% is small for all elements of q.
maxiter = 500;
iter = 0;
while any(abs(h(:)) > reltol*q(:))
    iter = iter + 1;
    if iter > maxiter
        % Too many iterations.  This should not happen.
        didnt = find(abs(h) > reltol*q); didnt = didnt(1);
        if numel(mu) == 1, mubad = mu; else mubad = mu(didnt); end
        if numel(lambda) == 1, lambdabad = b; else lambdabad = lambda(didnt); end
        if numel(p) == 1, pbad = p; else pbad = p(didnt); end
        warning('stats:addinvg:NoConvergence',...
                'INVGINV did not converge for mu = %g, lambda = %g, p = %g.',...
                mubad,lambdabad,pbad);
        break
    end

    h = (invgcdf(q,1,lambda0) - p) ./ max(invgpdf(q,1,lambda0),realmin);
    qnew = q - h;
    % Make sure that the current iterates stay positive.  When Newton's
    % Method suggests steps that lead to negative values, take a step
    % 9/10ths of the way to zero instead.
    ksmall = find(qnew <= 0);
    if ~isempty(ksmall)
        qnew(ksmall) = q(ksmall) / 10;
        h = q - qnew;
    end
    q = qnew;
end

% Add in the scale factor, and broadcast the values to the correct place if
% need be.
if allOK
    x = q .* mu;
else
    x(k) = q .* mu;
end


function r = invgrnd(mu, lambda, varargin)
%INVGRND Random arrays from the inverse Gaussian distribution.
mu(mu <= 0) = NaN;
lambda(lambda <= 0) = NaN;

[err, sizeOut] = statsizechk(2,mu,lambda,varargin{:});
if err > 0
    error('stats:invgrnd:InconsistentSizes','Size information is inconsistent.');
end

c = mu.*chi2rnd(1,sizeOut);
r = (mu./(2.*lambda)) .* (2.*lambda + c - sqrt(4.*lambda.*c + c.^2));
invert = (rand(sizeOut).*(mu+r) > mu);
r(invert) = mu.^2 ./ r(invert);


function [m,v] = invgstat(mu, lambda)
%INVGSTAT Mean and variance for the inverse Gaussian distribution.
mu(mu <= 0) = NaN;
lambda(lambda <= 0) = NaN;

m = mu;
v = mu.^3 ./ lambda;


function [nlogL,acov] = invglike(params,data,cens,freq)
%INVGLIKE Negative log-likelihood for the inverse Gaussian distribution.
if nargin < 4 || isempty(freq), freq = ones(size(data)); end
if nargin < 3 || isempty(cens), cens = zeros(size(data)); end

nlogL = invg_nloglf(params, data, cens, freq);
if nargout > 1
    acov = mlecov(params, data, 'nloglf',@invg_nloglf, 'cens',cens, 'freq',freq);
end


% ==== inverse Gaussian fitting functions ====

function [phat,pci] = invgfit(x,alpha,cens,freq,opts)
%INVGFIT Parameter estimates and confidence intervals for inverse Gaussian data.

if nargin < 2 || isempty(alpha), alpha = .05; end
if nargin < 3 || isempty(cens), cens = zeros(size(x)); end
if nargin < 4 || isempty(freq), freq = ones(size(x)); end
if nargin < 5, opts = []; end

if any(x <= 0)
    error('stats:invgfit:BadData','The data in X must be positive');
end

ncen = sum(freq.*cens);
if ncen == 0
    xbar = mean(x);
    phat = [xbar 1./mean(1./x - 1./xbar)];

else
    % MLEs of the uncensored data as starting point
    xunc = x(cens == 0); xbarunc = mean(xunc);
    start = [xbarunc 1./mean(1./xunc - 1./xbarunc)];

    % The default options include turning statsfminbx's display off.  This
    % function gives its own warning/error messages, and the caller can
    % turn display on to get the text output from statsfminbx if desired.
    options = statset(statset('invgfit'), opts);
    tolBnd = options.TolBnd;
    options = optimset(options);
    dfltOptions = struct('DerivativeCheck','off', 'HessMult',[], ...
        'HessPattern',ones(2,2), 'PrecondBandWidth',Inf, ...
        'TypicalX',ones(2,1), 'MaxPCGIter',1, 'TolPCG',0.1);

    % Maximize the log-likelihood with respect to mu and lambda.
    funfcn = {'fungrad' 'invgfit' @invg_nloglf [] []};
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
        warning('stats:invgfit:IterOrEvalLimit',wmsg);
    elseif (err < 0)
        error('stats:invgfit:NoSolution',...
              'Unable to reach a maximum likelihood solution.');
    end
end

% Compute CIs using a normal approximation for phat.
if nargout > 1
    acov = mlecov(phat, x, 'nloglf',@invg_nloglf, 'cens',cens, 'freq',freq);
    probs = [alpha/2; 1-alpha/2];
    se = sqrt(diag(acov))';
    pci = norminv([probs probs], [phat; phat], [se; se]);
end


function [nll,ngrad] = invg_nloglf(params, x, cens, freq)
%INVG_NLOGLF Objective function for inverse Gaussian maximum likelihood.

mu = params(1);
lambda = params(2);

L = .5.*log(lambda) - 1.5.*log(x) - lambda.*(x./mu-1).^2 ./ (2.*x);
ncen = sum(freq.*cens);
if ncen > 0
    cen = (cens == 1);
    xcen = x(cen);
    tmpsqrt = sqrt(lambda./xcen);
    tmpexp = exp(2.*lambda./mu);
    zcen = -(xcen./mu-1) .* tmpsqrt;
    wcen = -(xcen./mu+1) .* tmpsqrt;
    Phizcen = 0.5.*erfc(-zcen./sqrt(2));
    Phiwcen = 0.5.*erfc(-wcen./sqrt(2));
    Scen = Phizcen - tmpexp .* Phiwcen;
    L(cen) = log(Scen);
end
nll = -sum(freq .* L);

if nargout > 1
    dL1 = lambda.*(x-mu)./mu.^3;
    dL2 = 1./(2.*lambda) - (x./mu-1).^2 ./ (2.*x);
    if ncen > 0
        phizcen = exp(-0.5.*zcen.^2)./sqrt(2.*pi);
        phiwcen = exp(-0.5.*wcen.^2)./sqrt(2.*pi);
        dS1cen = (phizcen - tmpexp.*phiwcen).*(xcen./mu.^2).*tmpsqrt ...
                  + 2.*Phiwcen.*tmpexp.*lambda./mu.^2;
        dS2cen = 0.5.*(phizcen.*zcen - tmpexp.*phiwcen.*wcen)./lambda ...
                  - 2.*Phiwcen.*tmpexp./mu;
        dL1(cen) = dS1cen ./ Scen;
        dL2(cen) = dS2cen ./ Scen;
    end
    ngrad = -[sum(freq .* dL1) sum(freq .* dL2)];
end
