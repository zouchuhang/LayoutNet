function x = stdrinv(p, v, r)
%STDRINV Compute inverse c.d.f. for Studentized Range statistic
%   STDRINV(P,V,R) is the inverse cumulative distribution function for
%   the Studentized range statistic for R samples and V degrees of
%   freedom, evaluated at P.

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.3 $  $Date: 2002/02/04 18:52:50 $

% Based on Fortran program from statlib, http://lib.stat.cmu.edu
% Algorithm AS 190  Appl. Statist. (1983) Vol.32, No. 2
% Incorporates corrections from Appl. Statist. (1985) Vol.34 (1)

if (length(p)>1 | length(v)>1 | length(r)>1),
   error('STDRINV requires scalar arguments.'); % for now
end

[err,p,v,r] = distchck(3,p,v,r);
if (err > 0), error('Non-scalar arguments must match in size.'); end

% Handle illegal or trivial values first.
x = zeros(size(p));
if (length(x) == 0), return; end
ok = (v>0) & (v==round(v)) & (r>1) & (r==round(r) & (p<1));
x(~ok) = NaN;
ok = ok & (p>0);
v = v(ok);
p = p(ok);
r = r(ok);
if (length(v) == 0), return; end
xx = zeros(size(v));

% Define constants
jmax = 20;
pcut = 0.00001;
tiny = 0.000001;
upper = (p > .99);
if (upper)
   uppertail = 'u';
   p0 = 1-p;
else
   uppertail = 'l';
   p0 = p;
end

% Obtain initial values
q1 = qtrng0(p, v, r);
p1 = stdrcdf(q1, v, r, uppertail);
xx = q1;
if (abs(p1-p0) >= pcut*p0)
   if (p1 > p0), p2 = max(.75*p0, p0-.75*(p1-p0)); end
   if (p1 < p0), p2 = p0 + (p0 - p1) .* (1 - p0) ./ (1 - p1) * 0.75; end
   if (upper)
      q2 = qtrng0(1-p2, v, r);
   else
      q2 = qtrng0(p2, v, r);
   end

   % Refine approximation
   for j=2:jmax
      p2 = stdrcdf(q2, v, r, uppertail);
      e1 = p1 - p0;
      e2 = p2 - p0;
      d = e2 - e1;
      xx = (q1 + q2) / 2;
      if (abs(d) > tiny*p0)
         xx = (e2 .* q1 - e1 .* q2) ./ d;
      end
      if (abs(e1) >= abs(e2))
         q1 = q2;
         p1 = p2;
      end
      if (abs(p1 - p0) < pcut*p0), break; end
	   q2 = xx;
   end
end
   
x(ok) = xx;

% ---------------------------------
function x = qtrng0(p, v, r)
% Algorithm AS 190.2  Appl. Statist. (1983) Vol.32, No.2
% Calculates an initial quantile p for a studentized range
% distribution having v degrees of freedom and r samples
% for probability p, p.gt.0.80 .and. p.lt.0.995.

t=norminv(0.5 + 0.5 .* p);
if (v < 120), t = t + 0.25 * (t.^3 + t) ./ v; end
q = 0.8843 - 0.2368 .* t;
if (v < 120), q = q - (1.214./v) + (1.208.*t./v); end
x = t .* (q .* log(r-1) + 1.4142);
