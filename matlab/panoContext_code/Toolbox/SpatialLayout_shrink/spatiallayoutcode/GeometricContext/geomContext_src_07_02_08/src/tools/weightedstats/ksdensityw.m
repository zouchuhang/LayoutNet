function [f,x,u]=ksdensityw(y,w,varargin)
%KSDENSITY Compute density estimate
%   [F,XI]=KSDENSITY(X) computes a probability density estimate of the sample
%   in the vector X.  F is the vector of density values evaluated at the
%   points in XI.  The estimate is based on a normal kernel function, using a
%   window parameter (bandwidth) that is a function of the number of points
%   in X.  The density is evaluated at 100 equally-spaced points covering
%   the range of the data in X.
%
%   F=KSDENSITY(X,XI) specifies the vector XI of values where the density
%   estimate is to be evaluated.
%
%   [F,XI,U]=KSDENSITY(...) also returns the bandwidth of the kernel smoothing
%   window.
%
%   DWH: w are the weights for y (w = 1/n for ksdenstiy)
%
%   [...]=KSDENSITY(...,'PARAM1',val1,'PARAM2',val2,...) specifies parameter
%   name/value pairs to control the density estimation.  Valid parameters
%   are the following:
%
%      Parameter    Value
%      'kernel'     The type of kernel smoother to use, chosen from among
%                   'normal' (default), 'box', 'triangle', and
%                   'epanechinikov'.
%      'npoints'    The number of equally-spaced points in XI.
%      'width'      The bandwidth of the kernel smoothing window.  The default
%                   is optimal for estimating normal densities, but you
%                   may want to choose a smaller value to reveal features
%                   such as multiple modes.
%
%   In place of the kernel functions listed above, you can specify another
%   function by using @ (such as @normpdf) or quotes (such as 'normpdf').
%   The function must take a single argument that is an array of distances
%   between data values and places where the density is evaluated, and 
%   return an array of the same size containing corresponding values of
%   the kernel function.
%
%   Example:
%      x = [randn(30,1); 5+randn(30,1)];
%      [f,xi] = ksdensity(x);
%      plot(xi,f);
%   This example generates a mixture of two normal distributions, and
%   plots the estimated density.
%
%   See also HIST, @.

% Reference:
%   A.W. Bowman and A. Azzalini (1997), "Applied Smoothing
%      Techniques for Data Analysis," Oxford University Press.

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.7 $  $Date: 2002/03/21 20:36:29 $
   
% Get y vector and its dimensions
if (prod(size(y)) > length(y)), error('X must be a vector'); end
y = y(:);
y(isnan(y)) = [];
n = length(y);
ymin = min(y);
ymax = max(y);

% Maybe x was specified, or maybe not
if ~isempty(varargin)
   if ~ischar(varargin{1})
      x = varargin{1};
      varargin(1) = [];
   end
end

% Process additional name/value pair arguments
okargs = {'width' 'npoints' 'kernel'};
defaults = {[] [] 'normal'};
[emsg,u,m,kernel] = statgetargs(okargs, defaults, varargin{:});
error(emsg);

% Default window parameter is optimal for normal distribution
if (isempty(u)),   
   med = median(y);
   sig = median(abs(y-med)) / 0.6745;
   if sig<=0, sig = ymax-ymin; end
   if sig>0
      u = sig * (4/(3*n))^(1/5);
   else
      u = 1;
   end   
end

% Check other arguments or get defaults.
if ~exist('x','var')
   if isempty(m), m=100; end
   x = linspace(ymin-2*u, ymax+2*u, m);
elseif (prod(size(x)) > length(x))
   error('XI must be a vector');
end
xsize = size(x);
x = x(:);
m = length(x);

okkernels = {'normal' 'epanechinikov' 'box' 'triangle'};
if isempty(kernel)
   kernel = okkernels{1};
elseif ~(isa(kernel,'function_handle') | isa(kernel,'inline'))
   if ~ischar(kernel)
      error('Smoothing kernel must be a function.');
   end
   knum = strmatch(lower(kernel), okkernels);
   if (length(knum) == 1)
      kernel = okkernels{knum};
   end
end

blocksize = 1e6;
if n*m<=blocksize
   % Compute kernel density estimate in one operation
   z = (repmat(x',n,1)-repmat(y,1,m))/u;
   w2 = repmat(w, 1, m);
   f = sum(feval(kernel, z).*w2,1);
else
   % For large vectors, process blocks of elements as a group
   M = max(1,ceil(blocksize/n));
   mrem = rem(m,M);
   if mrem==0, mrem = min(M,m); end
   x = x';
   
   f = zeros(1,m);
   ii = 1:mrem;
   z = (repmat(x(ii),n,1)-repmat(y,1,mrem))/u;
   w2 = repmat(w, 1, mrem);
   f(ii) = sum(feval(kernel, z).*w2,1);
   z = zeros(n,M);
   w2 = zeros(n,M);
   for j=mrem+1:M:m
      ii = j:j+M-1;
      z(:) = (repmat(x(ii),n,1)-repmat(y,1,M))/u;
      w2(:) = repmat(w, 1, M);
      f(ii) = sum(feval(kernel, z).*w2,1);
   end
end

f = reshape(f./u, xsize);

% -----------------------------
% The following are functions that define smoothing kernels k(z).
% Each function takes a single input Z and returns the value of
% the smoothing kernel.  These sample kernels are designed to
% produce outputs that are somewhat comparable (differences due
% to shape rather than scale), so they are all probability
% density functions with unit variance.
%
% The density estimate has the form
%    f(x;k,h) = mean over i=1:n of k((x-y(i))/h) / h

function f = normal(z)
%NORMAL Normal density kernel.
%f = normpdf(z);
f = exp(-0.5 * z .^2) ./ sqrt(2*pi);

function f = epanechinikov(z)
%EPANECHINIKOV Epanechinikov's asymptotically optimal kernel.
a = sqrt(5);
z = max(-a, min(z,a));
f = .75 * (1 - .2*z.^2) / a;

function f = box(z)
%BOX    Box-shaped kernel
a = sqrt(3);
f = (abs(z)<=a) ./ (2 * a);

function f = triangle(z)
%TRIANGLE Triangular kernel.
a = sqrt(6);
z = abs(z);
f = (z<=a) .* (1 - z/a) / a; 