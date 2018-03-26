function [dy,y,d2y] = dgammainc(x,a,tail)
%DGAMMAINC Incomplete gamma function with derivatives.
%   DY = DGAMMAINC(X,A) computes the first derivative of the incomplete
%   gamma function, with respect to its second argument, for corresponding
%   elements of X and A.  X and A must be real and the same size (or either
%   can be a scalar).  A must also be non-negative.
%
%   [DY,Y] = DGAMMAINC(X,A) also returns the incomplete gamma function
%   itself.
%
%   [DY,Y,D2Y] = DGAMMAINC(X,A) also returns the second derivative of the
%   incomplete gamma function with respect to its second argument.
%
%   The incomplete gamma function is defined as:
%
%    gammainc(x,a) = 1 ./ gamma(a) .*
%       integral from 0 to x of t^(a-1).*exp(-t) dt
%
%   For any a >= 0, as x approaches infinity, y approaches 1 and both
%   derivatives approach 0.  For small x and a, y is approximately x^a, so
%   dgammainc(0,0) returns y == 1, dy == -Inf, and d2y == Inf.
%
%   [..] = GAMMAINC(X,A,TAIL) specifies the tail of the incomplete gamma
%   function when X is non-negative.  Choices are 'lower' (the default) and
%   'upper'.  The upper incomplete gamma function is defined as
%   1 - gammainc(x,a).
%
%   Warning: When x is negative, results can be inaccurate for abs(x) > a+1.
%
%   See also GAMMA, GAMMALN, GAMMAINC, PSI.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2004/01/24 09:36:09 $

if nargin < 3
    lower = true;
else
    switch tail
    case 'lower', lower = true;
    case 'upper', lower = false;
    otherwise, error('stats:dgammainc:InvalidTailArg', ...
                     'TAIL must be ''lower'' or ''upper''.');
    end
end

% x and a must be compatible for addition.
try
   y = x + a;
   y(:) = NaN;
catch
   error('stats:dgammainc:InputSizeMismatch', ...
         'X and A must be the same size, or scalars.')
end
dy = y;
if nargout > 2
    d2y = y;
end

if any(a(:) < 0)
   error('stats:dgammainc:NegativeArg', 'A must be non-negative.')
end

% If a is a vector, make sure x is too.
ascalar = isscalar(a);
if ~ascalar && isscalar(x)
   x = repmat(x,size(a));
end

% Upper limit for series and continued fraction.
amax = 2^20;

% Approximation for a > amax.  Accurate to about 5.e-5.
k = find(a > amax);
if ~isempty(k)
   if ascalar
      x = max(amax-1/3 + sqrt(amax/a).*(x-(a-1/3)),0);
      a = amax;
   else
      x(k) = max(amax-1/3 + sqrt(amax./a(k)).*(x(k)-(a(k)-1/3)),0);
      a(k) = amax;
   end
end

%
% Series expansion for lower incomplete gamma when x < a+1
%
k = find(x < a+1 & x ~= 0);
if ~isempty(k)
    xk = x(k);
    if ascalar, ak = a; else ak = a(k); end
    aplusn = ak;
    del = 1;
    ddel = 0;
    d2del = 0;
    sum = del;
    dsum = ddel;
    d2sum = d2del;
    while norm(del,'inf') >= 100*eps(norm(sum,'inf'))
        aplusn = aplusn + 1;
        del = del .* xk ./ aplusn;
        ddel = (ddel .* xk - del) ./ aplusn;
        d2del = (d2del .* xk - 2 .* ddel) ./ aplusn;
        sum = sum + del;
        dsum = dsum + ddel;
        d2sum = d2sum + d2del;
    end
    fac = exp(-xk + ak.*log(xk) - gammaln(ak+1));
    yk = fac.*sum;
    % For very small a, the series may overshoot very slightly.
    yk(xk > 0 & yk > 1) = 1;
    if lower, y(k) = yk; else y(k) = 1 - yk; end
    
    dlogfac = (log(xk) - psi(ak+1));
    dfac = fac .* dlogfac;
    dyk = dfac.*sum + fac.*dsum;
    if lower, dy(k) = dyk; else dy(k) = -dyk; end
    
    if nargout > 2
        d2fac = dfac.*dlogfac - fac.*psi(1,ak+1);
        d2yk = d2fac.*sum + 2.*dfac.*dsum + fac.*d2sum;
        if lower, d2y(k) = d2yk; else d2y(k) = -d2yk; end
    end
end

%
% Continued fraction for upper incomplete gamma when x >= a+1
%
k = find(x >= a+1); % & x ~= 0
if ~isempty(k)
    xk = x(k);
    if ascalar, ak = a; else ak = a(k); end
    n = 0;
    a0 = 0;
    a1 = ak;
    b0 = 1;
    b1 = xk;
    da0 = 0; db0 = 0; da1 = 1; db1 = 0;
    d2a0 = 0; d2b0 = 0; d2a1 = 0; d2b1 = 0;
    g = ak ./ xk;
    dg = 1 ./ xk;
    d2g = 0;
    d2gold = 1; % force one iteration, any nonzero value will do
    % Testing d2g is the more stringent than testing g or dg.  d2g may
    % be zero, so use a strict inequality
    while norm(d2g-d2gold,'inf') > 100*eps(norm(d2g,'inf'))
        rescale = 1 ./ b1; % keep terms from overflowing
        n = n + 1;
        nminusa = n - ak;
        d2a0 = (d2a1 + d2a0 .* nminusa - 2 .* da0) .* rescale;
        d2b0 = (d2b1 + d2b0 .* nminusa - 2 .* db0) .* rescale;
        da0 = (da1 + da0 .* nminusa - a0) .* rescale;
        db0 = (db1 + db0 .* nminusa - b0) .* rescale;
        a0 = (a1 + a0 .* nminusa) .* rescale;
        b0 = 1 + (b0 .* nminusa) .* rescale; % (b1 + b0 .* nminusa) .* rescale
        nrescale = n .* rescale;
        d2a1 = d2a0 .* xk + d2a1 .* nrescale;
        d2b1 = d2b0 .* xk + d2b1 .* nrescale;
        da1 = da0 .* xk + da1 .* nrescale;
        db1 = db0 .* xk + db1 .* nrescale;
        a1 = a0 .* xk + a1 .* nrescale;
        b1 = b0 .* xk + n; % b0 .* xk + b1 .* nrescale
        d2gold = d2g;
        g = a1 ./ b1;
        dg = (da1 - g.*db1) ./ b1;
        d2g = (d2a1 - dg.*db1 - g.*d2b1 - dg.*db1) ./ b1;
    end
    fac = exp(-xk + ak.*log(xk) - gammaln(ak+1));
    yk = fac.*g;
    if lower, y(k) = 1 - yk; else y(k) = yk; end
    
    dlogfac = (log(xk) - psi(ak+1));
    dfac = fac .* dlogfac;
    dyk = dfac.*g + fac.*dg;
    if lower, dy(k) = -dyk; else dy(k) = dyk; end
    
    if nargout > 2
        d2fac = dfac.*dlogfac - fac.*psi(1,ak+1);
        d2yk = d2fac.*g + 2.*dfac.*dg + fac.*d2g;
        if lower, d2y(k) = -d2yk; else d2y(k) = d2yk; end
    end
end

% Handle x == 0 separately to get it exactly correct.
kx0 = find(x == 0);
if ~isempty(kx0)
    if lower, y(kx0) = 0; else y(kx0) = 1; end
    dy(kx0) = 0;
    if nargout > 2
        d2y(kx0) = 0;
    end
end

% a == 0, x ~= 0 is already handled by the power series or continued
% fraction, now fill in dgammainc(0,0).  While we're at it, make
% gammainc(x,0) or 1-gammainc(x,0) exact for any x, not just x == 0.
ka0 = find(a == 0);
if ~isempty(ka0)
    if ascalar
        if lower
            y(:) = 1;
            dy(kx0) = -Inf;
            if nargout > 2, d2y(kx0) = Inf; end
        else
            y(:) = 0; 
            dy(kx0) = Inf;
            if nargout > 2, d2y(kx0) = -Inf; end
        end
    else
        ka0x0 = find(a == 0 & x == 0);
        if lower
            y(ka0) = 1;
            dy(ka0x0) = -Inf;
            if nargout > 2, d2y(ka0x0) = Inf; end
        else
            y(ka0) = 0;
            dy(ka0x0) = Inf;
            if nargout > 2, d2y(ka0x0) = -Inf; end
        end
    end
end
