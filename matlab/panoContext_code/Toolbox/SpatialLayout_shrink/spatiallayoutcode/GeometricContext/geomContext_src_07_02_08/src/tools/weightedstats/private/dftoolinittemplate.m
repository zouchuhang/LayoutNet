function s = dfittooldists
%DFITTOOLDISTS Initialize dfittool with custom distributions.
%
%    S=DFITTOOLDISTS is called during the initialization of DFITTOOL to get
%    any custom distributions you may want to define.  This function should
%    appear somewhere on your MATLAB path.  You can edit it to define
%    distributions that you want to be available for fitting in DFITTOOL.
%    S is a structure or an array of structures with fields as defined below.
%
%    You can load this after initialization by using the menu item
%    File -> Import Custom Distributions.  In that case you are not
%    restricted to use the name DFITTOOLDISTS.
%   
%    See also DFITTOOL.

%   Copyright 2001-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2004/01/24 09:35:58 $

% Create a structure to receive distribution information
s = struct;

% ---------------------------------------------------------
% ---- Remove the following return statement to define the 
% ---- Laplace distributon
% ---------------------------------------------------------
return

% Laplace (double exponential) distribution definition
j = 1;                                     % custom distribution #1
s(j).name = 'Laplace';                     % name for display
s(j).pnames = {'mu' 'sigma'};              % names of parameters
s(j).pdescription = {'location' 'scale'};  % descriptions of parameters
s(j).cdffunc = @laplcdf;                   % function to compute cdf
s(j).pdffunc = @laplpdf;                   % function to compute density
s(j).invfunc = @laplinv;                   % function to compute inverse cdf
s(j).fitfunc = @laplfit;                   % function to do fit
s(j).statfunc = @laplstat;                 % function to compute mean and var
s(j).islocscale = true;                    % location/scale distribution?


% ----------------------------------------------------------------
% ---- To enter your own distribution below, remove the following
% ---- return statement and modify the other statements as necessary
% ----------------------------------------------------------------
return

% Increment index
j = j+1;

% Enter a text string to display as distribution name
s(j).name = 'Enter Name Here';

% Enter the names of the parameters
s(j).pnames = {'p1', 'p2'};

% (optional) Enter a short description of each parameter (default empty)
%s(j).pdescription = {'location' 'scale'};

% (optional) Enter a vector indicating whether each parameter must have
% its value specified (true) or if it can be estimated (false, default)
%s(j).prequired = [false false];

% Enter function handles to compute the cdf, pdf, and inverse cdf
s(j).cdffunc = @yourcdf;
s(j).pdffunc = @yourpdf;
s(j).invfunc = @yourinv;

% (optional) Enter a function handle to compute the mean and variance
%s(j).statfunc = @yourstat;

% (optional) Enter a code name to use internally (default is lower case
% version of the distribution name)
%s(j).code = 'entercode';

% (optional) Is this a continuouse distribution (true, default) or is it
% on the integers only (false)
%s(j).iscontinuous = true;

% (optional) Is this distribution a location/scale family (default false)
%s(j).islocscale = false;

% (optional) Define a function that can fit this distribution (default none)
%s(j).fitfunc = @yourfit;

% (optional) Define a function that can compute the negative log
% likelihood (default none)
%s(j).likefunc = @yourlike;

% (optional) Define functions that can compute the cdf and inverse cdf
% on the log scale, for example @normcdf and @norminv can do this for
% the lognormal distribution (default none)
%s(j).logcdffunc = @yourlogcdf;
%s(j).loginvfunc = @yourloginv;

% (optional) Do the cdf and inverse cdf functions return confidence bounds
% as additional outputs? (default false)
%s(j).hasconfbounds = false;

% (optional) Does the fit function support censoring? (default false)
%s(j).censoring = false;

% (optional) Does the fit function return the fitted parameters as a
% single vector (true, default) or as separate scalars (false)
%s(j).paramvec = true;

% (optional) Enter a two-element vector defining the range over which
% this distribution gives positive probability, such as:
%    [-Inf Inf]    The entire real line (default)
%    [0 Inf]       Positive values only
%    [-1 1]        Values between -1 and 1 only
%s(j).support = [-Inf Inf];

% (optional) Enter a two-element vector specifying whether a data value
% can be exactly at the boundary (true) or must be strictly within the
% boundary (false, default for both end points)
%s(j).closedbound = [false false];

% (optional) Should the probability plot be on the log scale? (default false)
%s(j).uselogpp = false;


% ---------------------------------------------------
% ---- Define additonal distributions similarly below
% ---------------------------------------------------



% -----------------------------------------------------------------------
% ---- Include any local functions required for your custom distributions
% -----------------------------------------------------------------------

function p=laplfit(x,alpha)
%LAPLFIT Fit the Laplace distribution

% Fitting for this distribution is simple:
%    mu = median(x);
%    sigma = mean(abs(x-mu));
%    p = [mu sigma];

% For illustration, though, the following lines use a nonlinear
% fitting function that could also be used to fit other distributions

% Starting values are needed:
mu = mean(x);
sigma = std(x);

% Fit the Laplace distribution
p = mle(x,'pdf',@laplpdf,'cdf',@laplcdf,'start',[mu, sigma]);


function f=laplpdf(x,mu,sigma)
%LAPLDF Laplace distribution probability density function
if nargin<2, mu=0; end
if nargin<3, sigma=1; end
z = (x-mu)./sigma;
f = exp(-abs(z))/(2*sigma);


function f=laplcdf(x,mu,sigma)
%LAPLCDF Laplace distribution cumulative distribution function
if nargin<2, mu=0; end
if nargin<3, sigma=1; end
f = zeros(size(x));
z = (x-mu)./sigma;
t = (z<=0);
f(t) = exp(z(t))/2;
f(~t) = 1 - exp(-z(~t))/2;


function z=laplinv(p,mu,sigma)
%LAPLINV Laplace distribution inverse cumulative distribution function
if nargin<2, mu=0; end
if nargin<3, sigma=1; end
z = zeros(size(p));
t = (p<=.5);
z(t) = log(2*p(t));
z(~t) = -log(1-2*(p(~t)-.5));
z = mu + sigma*z;


function [m,v]=laplstat(mu,sigma)
%LAPLSTAT Laplace distribution statistics
m = mu;
v = 2 * sigma.^2;
          
