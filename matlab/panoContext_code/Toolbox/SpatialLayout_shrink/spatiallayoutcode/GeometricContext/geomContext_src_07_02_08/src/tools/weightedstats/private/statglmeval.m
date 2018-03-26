function retval=statglmeval(action,fn,varargin)
%STATGLMEVAL  Evaluate or test link function in proper environment
%   STATGLMEVAL('eval',FN,ARGS,...) evaluates the function FN in an
%   environment in which certain functions such as LOGIT and
%   D_LOGIT are defined.  This allows the function FN to be either
%   a user-defined function or one of the pre-defined functions
%   provided in the GLMFIT function, without contaminating the name
%   space with those pre-defined functions.
%
%   STATGLMEVAL('testlink',FN) test for the existence of FN as a
%   function handle, inline function, or text string containing the
%   name of a function in an M-file.

%   Author:  Tom Lane, 3-6-2000
%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.3 $  $Date: 2002/02/04 19:25:46 $

switch(action)
 case 'eval'
   retval = feval(fn,varargin{:});

 case 'testlink'
   c = class(fn);
   fnclass = class(@logit);
   if (isequal(c,fnclass) | isequal(c,'inline') | ...
       (isequal(c,'char') & ~isempty(which(fn))))
      retval = 1;
   else
      retval = 0;
   end
end
   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following functions are related to the link function that
% links the distribution parameter mu with the linear combination
% eta of predictor variables.  For example, the logit link defines
%
%        eta = log(mu/(1-mu))

%%%% functions for identity link
function a=identity(b)
a=b;

function a=d_identity(b)
a = ones(size(b));

function b=i_identity(a)
b=a;

%%%% functions for logit link
function a=logit(p)
a = log(p ./ (1-p));

function a=d_logit(p)
a = 1 ./ max(eps, (p .* (1-p)));

function p=i_logit(a)
p = 1 ./ (1 + exp(-a));


%%%% functions for probit link
function a=probit(p)
a = norminv(p);

function a=d_probit(p)
a = 1 ./ max(eps, normpdf(norminv(p)));

function p=i_probit(a)
p = normcdf(a);


%%%% functions for complementary log log link
function a=comploglog(p)
a = log(-log(1-max(eps,p)));

function a=d_comploglog(p)
a = 1 ./ -(max(eps,1-p) .* log(1-max(eps,p)));

function p=i_comploglog(a)
p = 1 - exp(-exp(a));


%%%% functions for log log link
function a=logloglink(p)
a = log(-log(max(eps,p)));

function a=d_logloglink(p)
a = 1 ./ (max(eps, p) .* log(max(eps,p)));

function p=i_logloglink(a)
p = exp(-exp(a));


%%%% functions for log link
function a=d_log(b);
a = 1 ./ max(eps,b);

function b=i_log(a);
b = exp(a);

%%%% functions for reciprocal link
function a=reciprocal(b)
a = 1 ./ max(eps, b);

function a=d_reciprocal(b);
a = -1 ./ max(eps,b).^2;

function b=i_reciprocal(a);
b = 1 ./ max(eps, a);


%%%% functions for power link
function a=power(b,p)
if (p==0)
   a = log(max(eps,b));
else
   a = max(eps,b) .^ p;
end

function a=d_power(b,p);
if (p==0)
   a = 1 ./ max(eps,b);
else
   a = p * max(eps,b).^(p-1);
end

function b=i_power(a,p);
if (p==0)
   b = exp(a);
else
   b = max(eps,a) .^ (1/p);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following functions define the variance 

function a=normalvariance(b)
a = ones(size(b));

function a=poissonvariance(b)
a = b;

function a=binomialvariance(p,N)
a = p .* (1-p) ./ N;

function a=gammavariance(b)
a = b.^2;

function a=inversegaussianvariance(b)
a = b.^3;
   