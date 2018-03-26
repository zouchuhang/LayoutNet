function xout = stdrcdf(q, v, r, upper)
%STDRCDF Compute c.d.f. for Studentized Range statistic
%   F = STDRCDF(Q,V,R) is the cumulative distribution function for the
%   Studentized range statistic for R samples and V degrees of
%   freedom, evaluated at Q.
%
%   G = STDRCDF(Q,V,R,'upper') is the upper tail probability,
%   G=1-F.  This version computes the upper tail probability
%   directly (not by subtracting it from 1), and is likely to be
%   more accurate if Q is large and therefore F is close to 1.

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.3 $  $Date: 2002/02/04 19:25:50 $

% Based on Fortran program from statlib, http://lib.stat.cmu.edu
% Algorithm AS 190  Appl. Statist. (1983) Vol.32, No. 2
% Incorporates corrections from Appl. Statist. (1985) Vol.34 (1)
% Vectorized and simplified for MATLAB.  Added 'upper' option.

if (length(q)>1 | length(v)>1 | length(r)>1),
   error('STDRCDF requires scalar arguments.'); % for now
end
[err,q,v,r] = distchck(3,q,v,r);
if (err > 0), error('Non-scalar arguments must match in size.'); end
uppertail = 0;
if (nargin>3)
   if ~(  isequal(upper,'u') | isequal(upper,'upper') ...
        | isequal(upper,'l') | isequal(upper,'lower'))
      error('Fourth argument must be ''upper'' or ''lower''.');
   end
   uppertail = isequal(upper,'u') | isequal(upper,'upper');
end

% Accuracy can be increased by use of a finer grid.  Increase
% jmax, kmax and 1/step proportionally.
jmax = 15;          % controls maximum number of steps
kmax = 15;          % controls maximum number of steps
step = 0.45;        % node spacing
vmax = 120;         % max d.f. for integration over chi-square

% Handle illegal or trivial values first.
xout = zeros(size(q));
if (length(xout) == 0), return; end   
ok = (v>0) & (v==round(v)) & (r>1) & (r==round(r));
xout(~ok) = NaN;
ok = ok & (q > 0);
v = v(ok);
q = q(ok);
r = r(ok);
if (length(v) == 0), return; end
xx = zeros(size(v));

% Compute constants, locate midpoint, adjust steps.
g = step ./ (r .^ 0.2);
if (v > vmax)
   c = log(r .* g ./ sqrt(2*pi));
else
   h = step ./ sqrt(v);
   v2 = v * 0.5;
   c = sqrt(2/pi) * exp(-v2) .* (v2.^v2) ./ gamma(v2);
   c = log(c .* r .* g .* h);

   j=(-jmax:jmax)';
   hj = h * j;
   ehj = exp(hj);
   qw = q .* ehj;
   vw = v .* (hj + 0.5 * (1 - ehj .^2));
   C = ones(1,2*kmax+1);         % index to duplicate columns
   R = ones(1,2*jmax+1);         % index to duplicate rows
end

% Compute integral by summing the integrand over a
% two-dimensional grid centered approximately near its maximum.
gk = (0.5 * log(r)) + g * (-kmax:kmax);
w0 = c - 0.5 * gk .^ 2;
pz = normcdf(-gk);
if (~uppertail)
   % For regular cdf, use integrand as in AS 190.
   if (v > vmax)
      % don't integrate over chi-square
      x = normcdf(q - gk) - pz;
      xx = sum(exp(w0) .* (x .^ (r-1)));
   else
      % integrate over chi-square
      x = normcdf(qw(:,C) - gk(R,:)) - pz(R,:);
      xx = sum(sum(exp(w0(R,:) + vw(:,C)) .* (x .^ (r-1))));
   end
else
   % To compute the upper tail probability, we need an integrand that
   % contains the normal probability of a region consisting of a
   % hyper-quadrant minus a rectangular region at the origin of the
   % hyperquadrant.
   if (v > vmax)           % for large d.f., don't integrate over chi-square
      xhq   = (1 - pz) .^ (r-1);
      xrect = (normcdf(q - gk) - pz) .^ (r-1);
      xx = sum(exp(w0) .* (xhq - xrect));
   else                    % for typical cases, integrate over chi-square
      xhq   = (1 - pz) .^ (r-1);
      xrect = (normcdf(qw(:,C) - gk(R,:)) - pz(R,:)) .^ (r-1);
      xx = sum(sum(exp(w0(R,:) + vw(:,C)) .* (xhq(R,:) - xrect)));
   end
end

xout(ok) = xx;