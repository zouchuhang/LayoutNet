function s = addrice(s)
%ADDRICE Add the Rician distribution.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.10 $  $Date: 2004/02/01 22:10:34 $

j = length(s) + 1;
s(j).name = 'Rician';
s(j).code = 'rician';
s(j).pnames = {'s' 'sigma'};
s(j).pdescription = {'noncentrality' 'scale'};
s(j).prequired = [false false];
s(j).fitfunc = @ricefit;
s(j).likefunc = @ricelike;
s(j).cdffunc = @ricecdf;
s(j).pdffunc = @ricepdf;
s(j).invfunc = @riceinv;
s(j).statfunc = @ricestat;
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


% ==== Rician distribution functions ====

% these distribution functions do not yet handle arrays of parameters

function y = ricepdf(x,s,sigma)
%RICEPDF Rician probability density function (pdf).
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

x(x<0) = 0;
sigsq = sigma.^2;
expon = (x.^2 + s.^2)./(2.*sigsq);
y = (x./sigsq) .* exp(-expon) .* besseli(0, x.*s./sigsq);
y(expon > (log(realmax(class(x)))-1)) = 0; % fix up 0*Inf


function p = ricecdf(x,s,sigma)
%RICECDF Rician cumulative distribution function (cdf).
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

x(x<0) = 0;
p = ncx2cdf((x./sigma).^2, 2, (s./sigma).^2);


function x = riceinv(p,s,sigma)
%RICEINV Inverse of the Rician cumulative distribution function (cdf).
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

x = sigma .* sqrt(ncx2inv(p, 2, (s./sigma).^2));


function r = ricernd(s,sigma,varargin)
%RICERND Random arrays from the Rician distribution.
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

[err, sizeOut] = statsizechk(2,s,sigma,varargin{:});
if err > 0
    error('stats:ricernd:InconsistentSizes','Size information is inconsistent.');
end

r = sigma .* sqrt(ncx2rnd(2, (s./sigma).^2, sizeOut));


function [m,v] = ricestat(s,sigma)
%RICESTAT Mean and variance for the Rician distribution.
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

t = .5 .* (s./sigma).^2;
m = sigma.*sqrt(.5.*pi).*exp(-.5.*t) .* ((1+t).*besseli(0,.5.*t) + t.*besseli(1,.5.*t));
v = 2.*sigma.^2 + s.^2 - m.^2;


function [nlogL,acov] = ricelike(params,data,cens,freq)
%RICELIKE Negative log-likelihood for the Rician distribution.
if nargin < 4 || isempty(freq), freq = ones(size(data)); end
if nargin < 3 || isempty(cens), cens = zeros(size(data)); end

nlogL = rice_nloglf(params, data, cens, freq);
if nargout > 1
    acov = mlecov(params, data, 'nloglf',@rice_nloglf, 'cens',cens, 'freq',freq);
end


% ==== Rician fitting functions ====

function [phat,pci] = ricefit(x,alpha,cens,freq,opts)
%NAKAFIT Parameter estimates and confidence intervals for Rician data.

if nargin < 2 || isempty(alpha), alpha = .05; end
if nargin < 3 || isempty(cens), cens = zeros(size(x)); end
if nargin < 4 || isempty(freq), freq = ones(size(x)); end
if nargin < 5, opts = []; end

if any(x <= 0)
    error('stats:ricefit:BadData','The data in X must be positive');
end

% Moment estimators of the uncensored data as starting point
% E[x.^2] = s.^2 + 2.*sigma.^2
% E[x.^4] = s.^4 + 8.*s.^2.*sigma.^2 + 8.*sigma.^4
xsqunc = x(cens == 0).^2;
meanxsq = mean(xsqunc); meanx4th = mean(xsqunc);
if meanxsq.^2 < meanx4th && meanx4th < 2.*meanxsq.^2
    s4th = 2.*meanxsq.^2 - meanx4th;
    ssq = sqrt(s4th);
    sigsq = .5.*(meanxsq - ssq);
    start = [sqrt(ssq) sqrt(sigsq)];
else
    start = cast([1 1],class(x));
end

% The default options include turning fminsearch's display off.  This
% function gives its own warning/error messages, and the caller can turn
% display on to get the text output from fminsearch if desired.
options = statset(statset('ricefit'), opts);
tolBnd = options.TolBnd;
options = optimset(options);

% Maximize the log-likelihood with respect to s and sigma.
[phat,nll,err,output] = ...
    fminsearch(@rice_nloglf, start, options, x, cens, freq, tolBnd);
if (err == 0)
    % fminsearch may print its own output text; in any case give something
    % more statistical here, controllable via warning IDs.
    if output.funcCount >= options.MaxFunEvals
        wmsg = 'Maximum likelihood estimation did not converge.  Function evaluation limit exceeded.';
    else
        wmsg = 'Maximum likelihood estimation did not converge.  Iteration limit exceeded.';
    end
    warning('stats:ricefit:IterOrEvalLimit',wmsg);
elseif (err < 0)
    error('stats:ricefit:NoSolution',...
          'Unable to reach a maximum likelihood solution.');
end

% Compute CIs using a normal approximation for phat.
if nargout > 1
    acov = mlecov(phat, x, 'nloglf',@rice_nloglf, 'cens',cens, 'freq',freq);
    probs = [alpha/2; 1-alpha/2];
    se = sqrt(diag(acov))';
    pci = norminv([probs probs], [phat; phat], [se; se]);
end


function nll = rice_nloglf(parms, x, cens, freq, tolBnd)
%RICE_NLOGLF Objective function for Rician maximum likelihood.
s = parms(1);
sigma = parms(2);
sigsq = sigma.^2;

% Restrict sigma to the open interval (0, Inf).
if nargin > 4
    if s < tolBnd || sigma < tolBnd
        nll = Inf;
        return
    end
end

bess0 = besseli(0, x.*s./sigsq);
rsq = (x.^2 + s.^2)./(2.*sigsq);
L = -rsq + log(bess0) + log(x./sigsq);
ncen = sum(freq.*cens);
if ncen > 0
    cen = (cens == 1);
    xcen = x(cen);
    L(cen) = log(marcumq(s./sigma,xcen./sigma));
end
nll = -sum(freq .* L);

% Don't have derivatives of the Marcum's Q, so can't compute an analytic
% gradient with censoring.
%
% if nargout > 1
%     dlogbess0 = besseli(1, x.*s./sigsq) ./ bess0;
%     dL1 = (-s + dlogbess0.*x) ./ sigsq;
%     dL2 = (rsq - 1 - dlogbess0.*x.*s./sigsq) ./ sigma;
%     if ncen > 0
% %         dL1(cen) = ;
% %         dL2(cen) = ;
%     end
%     ngrad = -[sum(freq .* dL1) sum(freq .* dL2)];
% end


function Q = marcumq(a,b)
% Q = MARCUMQ(A,B) returns Marcum's "Q" function.

if isa(a,'single') || isa(b,'single')
   Q = repmat(single(NaN), size(a));
else
   Q = repmat(NaN, size(a));
end
Q(a~=Inf & b==0) = 1;
Q(a~=Inf & b==Inf) = 0;
Q(a==Inf & b~=Inf) = 1;
z = (isnan(Q) & a==0 & b~=Inf);
if (any(z))
   Q(z) = exp((-b(z).^2)./2);
end

z = isnan(Q) & ~isnan(a) & ~isnan(b);
if (any(z(:)))
%    aa = (a(z).^2)./2;
   aa = (a.^2)./2;
   bb = (b(z).^2)./2;

   d = exp(-aa);
   h = d;
   f = bb.*exp(-bb);
   k = 1;
   delta = f .* h;
   sum = delta;
   j = (delta > sum.*eps(class(delta)));
   while any(j)
      d = aa.*d./k;
      h = h + d;
      f = bb.*f./(k+1);
      delta = f .* h;
      sum(j) = sum(j) + delta(j);
      j = (delta > sum.*eps(class(delta)));
      k = k + 1;
   end
   Q(z) = 1 - sum;
end
