function[x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN]=statsfminbx(funfcn,x,l,u,verb,options,defaultopt,...
   computeLambda,initialf,initialGRAD,initialHESS,Hstr,varargin)
%SFMINBX Nonlinear minimization with box constraints.
%
% Locate a local minimizer to
%
%               min { f(x) :  l <= x <= u}.
%
%	where f(x) maps n-vectors to scalars.

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.5 $  $Date: 2002/06/26 18:22:31 $

%   Initialization
xcurr = x(:);  % x has "the" shape; xcurr is a vector
n = length(xcurr);
it= 1; totls = 0;

header = sprintf(['\n                                Norm of      First-order \n',...
      ' Iteration        f(x)          step          optimality   CG-iterations']);
formatstr = ' %5.0f      %13.6g  %13.6g   %12.3g     %7.0f';
if n == 0,
   error('n must be positive'),
end
if isempty(l),
   l = -inf*ones(n,1);
end,
if isempty(u),
   u = inf*ones(n,1);
end
arg = (u >= 1e10); arg2 = (l <= -1e10);
u(arg) = inf*ones(length(arg(arg>0)),1);
l(arg2) = -inf*ones(length(arg2(arg2>0)),1);
if any(u == l)
   errmsg=sprintf('%s\n%s',...
      'Equal upper and lower bounds not permitted in this large-scale method.',...
      'Use equality constraints and the medium-scale method instead.');
   error(errmsg)
elseif min(u-l) <= 0
   error('Inconsistent bounds.')
end
if min(min(u-xcurr),min(xcurr-l)) < 0, xcurr = startx(u,l); end

% get options out
typx = optimget(options,'TypicalX',defaultopt,'fast') ;
% In case the defaults were gathered from calling: optimset('quadprog'):
numberOfVariables = n;
if ischar(typx)
   if isequal(lower(typx),'ones(numberofvariables,1)')
      typx = ones(numberOfVariables,1);
   else
      error('Option ''TypicalX'' must be an integer value if not the default.')
   end
end

% Will be user-settable later:
pcmtx = optimget(options,'Preconditioner','hprecon') ; % not a default yet

mtxmpy = optimget(options,'HessMult',defaultopt,'fast') ;
if isempty(mtxmpy)
    mtxmpy = @hmult; % to detect name clash with user hmult.m, need this
end

active_tol = optimget(options,'ActiveConstrTol',sqrt(eps)); % not a default yet, so use slow optimget
pcflags = optimget(options,'PrecondBandWidth',defaultopt,'fast') ;
tol2 = optimget(options,'TolX',defaultopt,'fast') ;
tol1 = optimget(options,'TolFun',defaultopt,'fast') ;
tol = tol1;
maxiter = optimget(options,'MaxIter',defaultopt,'fast') ;
maxfunevals = optimget(options,'MaxFunEvals',defaultopt,'fast') ;
pcgtol = optimget(options,'TolPCG',defaultopt,'fast') ;  % pcgtol = .1;
kmax = optimget(options,'MaxPCGIter', defaultopt,'fast') ;
if ischar(kmax)
   if isequal(lower(kmax),'max(1,floor(numberofvariables/2))')
      kmax = max(1,floor(numberOfVariables/2));
   else
      error('Option ''MaxPCGIter'' must be an integer value if not the default.')
   end
end
if ischar(maxfunevals)
   if isequal(lower(maxfunevals),'100*numberofvariables')
      maxfunevals = 100*numberOfVariables;
   else
      error('Option ''MaxFunEvals'' must be an integer value if not the default.')
   end
end
maxcount = min(maxiter, maxfunevals); % numfunevals = iterations, so just take minimum

dnewt = []; gopt = [];
ex = 0; posdef = 1; npcg = 0;

%tol1 = tol; tol2 = sqrt(tol1)/10;
if strcmp(optimget(options,'DerivativeCheck',defaultopt,'fast'),'on')
   warnstr = sprintf('%s\n%s\n', ...
      'Trust region algorithm does not currently check user-supplied gradients,', ...
      '  ignoring OPTIONS.DerivativeCheck.');
   warning(warnstr);
end

vpos(1,1) = 1; vpcg(1,1) = 0; nbnds = 1;
pcgit = 0; delta = 10;nrmsx = 1; ratio = 0; degen = inf;
if (all(u == inf) & all(l == -inf)) nbnds = 0; end
DS = speye(n);   v = zeros(n,1); dv = ones(n,1); del = 10*eps;
oval = inf;  g = zeros(n,1); newgrad = g; Z = [];

% Make x conform to the user's input x
x(:) = xcurr;
%   Evaluate f,g,  and H
if ~isempty(Hstr)  % use sparse finite differencing
   %[val,g] = feval(fname,x);
   switch funfcn{1}
   case 'fun'
      error('should not reach this')
   case 'fungrad'
      %[val,g(:)] = feval(funfcn{3},x,varargin{:});
      val = initialf; g(:) = initialGRAD;
   case 'fun_then_grad'
      % val = feval(funfcn{3},x,varargin{:});
      % g(:) = feval(funfcn{4},x,varargin{:});
      val = initialf; g(:) = initialGRAD;
   otherwise
      if isequal(funfcn{2},'fmincon')
         error('Undefined calltype in FMINCON');
      else
         error('Undefined calltype in FMINUNC');
      end
   end

   %      Determine coloring/grouping for sparse finite-differencing
   p = colmmd(Hstr)'; p = (n+1)*ones(n,1)-p; group = color(Hstr,p);
   % pass in the user shaped x
   H = sfd(x,g,Hstr,group,[],funfcn,varargin{:});
   %
else % user-supplied computation of H or dnewt
   % [val,g,H] = feval(fname,x);
   switch funfcn{1}
   case 'fungradhess'
     % [val,g(:),H] = feval(funfcn{3},x,varargin{:});
      val = initialf; g(:) = initialGRAD; H = initialHESS;
   case 'fun_then_grad_then_hess'
      % val = feval(funfcn{3},x,varargin{:});
      % g(:) = feval(funfcn{4},x,varargin{:});
      % H = feval(funfcn{5},x,varargin{:});
      val = initialf; g(:) = initialGRAD; H = initialHESS;
   otherwise
      if isequal(funfcn{2},'fmincon')
         error('Undefined calltype in FMINCON');
      else
         error('Undefined calltype in FMINUNC');
      end
   end
end

delbnd = max(100*norm(xcurr),1);
[nn,pp] = size(g);

%   Extract the Newton direction?
if pp == 2, dnewt = g(1:n,2); end
if verb > 1
   disp(header)
end

%   MAIN LOOP: GENERATE FEAS. SEQ.  xcurr(it) S.T. f(x(it)) IS DECREASING.
while ~ex
   if ~isfinite(val) | any(~isfinite(g))
      errmsg= sprintf('%s%s%s',funfcn{2},' cannot continue: ',...
         'user function is returning Inf or NaN values.');
      error(errmsg)
   end

   %     Update
   [v,dv] = definev(g(:,1),xcurr,l,u);
   gopt = v.*g(:,1); gnrm = norm(gopt,inf);
   vgnrm(it,1)=gnrm;
   r = abs(min(u-xcurr,xcurr-l)); degen = min(r + abs(g(:,1)));
   vdeg(it,1) = min(degen,1); bndfeas = min(min(xcurr-l,u-xcurr));
   if ((u == inf*ones(n,1)) & (l == -inf*ones(n,1))) degen = -1; end

   % Display
   if verb > 1
      currOutput = sprintf(formatstr,it,val,nrmsx,gnrm,pcgit);
      disp(currOutput);
   end

   %     TEST FOR CONVERGENCE
   diff = abs(oval-val);
   oval = val;
   if (nrmsx < .9*delta)&(ratio > .25)&(diff < tol1*(1+abs(oval)))
      ex = 1;
      if verb > .5
         disp('Optimization terminated successfully:')
         disp(' Relative function value changing by less than OPTIONS.TolFun');
      end

   elseif (it > 1) & (nrmsx < tol2)
      ex = 2;
      if verb > .5
         disp('Optimization terminated successfully:')
         disp(' Norm of the current step is less than OPTIONS.TolX');
      end

   elseif ((gnrm < tol1) & (posdef ==1) )
      ex = 3;
      if verb > .5
         disp('Optimization terminated successfully:')
         disp(' First-order optimality less than OPTIONS.TolFun, and no negative/zero curvature detected');
      end
   end

   %     Step computation
   if ~ex

      %       Determine trust region correction
      dd = abs(v); D = sparse(1:n,1:n,full(sqrt(dd)));
      sx = zeros(n,1); theta = max(.95,1-gnrm);
      oposdef = posdef;
      [sx,snod,qp,posdef,pcgit,Z] = trdog(xcurr, g(:,1),H,D,delta,dv,...
         mtxmpy,pcmtx,pcflags,pcgtol,kmax,theta,l,u,Z,dnewt,'hessprecon',varargin{:});
      if isempty(posdef), posdef = oposdef; end
      nrmsx=norm(snod); npcg=npcg + pcgit;
      newx=xcurr + sx; vpcg(it+1,1)=pcgit;

      %       Perturb?
      [pert,newx] = perturb(newx,l,u);
      vpos(it+1,1) = posdef;

      % Make newx conform to user's input x
      x(:) = newx;
      % Evaluate f, g, and H
      if ~isempty(Hstr) % use sparse finite differencing
         %[newval,newgrad] = feval(fname,x);
         switch funfcn{1}
         case 'fun'
            error('should not reach this')
         case 'fungrad'
            [newval,newgrad(:)] = feval(funfcn{3},x,varargin{:});
         case 'fun_then_grad'
            newval = feval(funfcn{3},x,varargin{:});
            newgrad(:) = feval(funfcn{4},x,varargin{:});
         otherwise
            error('Undefined calltype in FMINUNC');
         end
         newH = sfd(x,newgrad,Hstr,group,[],funfcn,varargin{:});

      else % user-supplied computation of H or dnewt
         %[newval,newgrad,newH] = feval(fname,x);
         switch funfcn{1}
         case 'fungradhess'
            [newval,newgrad(:),newH] = feval(funfcn{3},x,varargin{:});
         case 'fun_then_grad_then_hess'
            newval = feval(funfcn{3},x,varargin{:});
            newgrad(:) = feval(funfcn{4},x,varargin{:});
            newH = feval(funfcn{5},x,varargin{:});
         otherwise
            error('Undefined calltype in FMINUNC');
         end

      end
      [nn,pp] = size(newgrad);
      aug = .5*snod'*((dv.*abs(newgrad(:,1))).*snod);
      ratio = (newval + aug -val)/qp; vratio(it,1) = ratio;

      if (ratio >= .75) & (nrmsx >= .9*delta)
         delta = min(delbnd,2*delta);
      elseif ratio <= .25
         delta = min(nrmsx/4,delta/4);
      end
      if newval == inf
         delta = min(nrmsx/20,delta/20);
      end

      %       Update
      if newval < val
         xold = xcurr; xcurr=newx; val = newval; g= newgrad; H = newH;
         Z = [];

         %          Extract the Newton direction?
         if pp == 2, dnewt = newgrad(1:n,2); end
      end
      it = it+1; vval(it,1) = val;
   end
   if it > maxcount,
      ex=4;
      it = it-1;
      if verb > 0
         if it > maxiter
            disp('Maximum number of iterations exceeded;')
            disp('   increase options.MaxIter')
         elseif it > maxfunevals
            disp('Maximum number of function evaluations exceeded;')
            disp('   increase options.MaxFunEvals')
         end
      end
   end
end % while

HESSIAN = H;
GRAD = g;
FVAL = val;
LAMBDA = [];
if ex==4
   EXITFLAG = 0;
elseif ex==10
   EXITFLAG = -1;
else
   EXITFLAG = 1;
end
OUTPUT.iterations = it;
OUTPUT.funcCount = it;
OUTPUT.cgiterations = npcg;
OUTPUT.firstorderopt = gnrm;
OUTPUT.algorithm = 'large-scale: trust-region reflective Newton';
x(:) = xcurr;
if computeLambda
   g = full(g);

   LAMBDA.lower = zeros(length(l),1);
   LAMBDA.upper = zeros(length(u),1);
   argl = logical(abs(xcurr-l) < active_tol);
   argu = logical(abs(xcurr-u) < active_tol);

   LAMBDA.lower(argl) = (g(argl));
   LAMBDA.upper(argu) = -(g(argu));
   LAMBDA.ineqlin = []; LAMBDA.eqlin = []; LAMBDA.ineqnonlin=[]; LAMBDA.eqnonlin=[];
else
   LAMBDA = [];
end


%===== definev.m =================================================


function [v,dv]= definev(g,x,l,u);
%DEFINEV Scaling vector and derivative
%
%	[v,dv]= DEFINEV(g,x,l,u) returns v, distances to the
%   bounds corresponding to the sign of the gradient g, where
%   l is the vector of lower bounds, u is the vector of upper
%   bounds. Vector dv is 0-1 sign vector (See ?? for more detail.)
%

n = length(x);
v = zeros(n,1);
dv=zeros(n,1);
arg1 = (g < 0)  & (u <  inf );
arg2 = (g >= 0) & (l > -inf);
arg3 = (g < 0)  & (u == inf);
arg4 = (g >= 0) & (l == -inf);
v(arg1)  = (x(arg1) - u(arg1));
dv(arg1) = 1;
v(arg2)  = (x(arg2) - l(arg2));
dv(arg2) = 1;
v(arg3)  = -1;
dv(arg3) = 0;
v(arg4)  = 1;
dv(arg4) = 0;


%===== trdog.m ===================================================


function[s,snod,qpval,posdef,pcgit,Z] = trdog(x,g,H,D,delta,dv,...
   mtxmpy,pcmtx,pcoptions,tol,kmax,theta,l,u,Z,dnewt,preconflag,varargin);
%TRDOG Reflected (2-D) trust region trial step (box constraints)
%
% [s,snod,qpval,posdef,pcgit,Z] = TRDOG(x,g,H,D,delta,dv,...
%                 mtxmpy,pcmtx,pcoptions,tol,theta,l,u,Z,dnewt,preconflag);
%
%   Determine the trial step `s', an approx. trust region solution.
%   `s' is chosen as the best of 3 steps: the scaled gradient
%   (truncated to  maintain strict feasibility),
%   a 2-D trust region solution (truncated to remain strictly feas.),
%   and the reflection of the 2-D trust region solution,
%   (truncated to remain strictly feasible).
%
%   The 2-D subspace (defining the trust region problem) is defined
%   by the scaled gradient direction and a CG process (returning
%   either an approximate Newton step of a direction of negative curvature.
%   Driver functions are: SNLS, SFMINBX
%   SNLS actually calls TRDOG with the Jacobian matrix (and a special
%   Jacobian-matrix multiply function in MTXMPY).

% Initialization
n = length(g);
pcgit = 0;
grad = D*g;
DM = D;
DG = sparse(1:n,1:n,full(abs(g).*dv));
posdef = 1;
pcgit = 0;
tol2 = sqrt(eps);
v1 = dnewt;
qpval1 = inf;
qpval2 = inf;
qpval3 = inf;

% DETERMINE A 2-DIMENSIONAL SUBSPACE
if isempty(Z)
   if isempty(v1)
      switch preconflag
      case 'hessprecon'
         % preconditioner based on H, no matter what it is
         [R,permR] = feval(pcmtx,H,pcoptions,DM,DG,varargin{:});
      case 'jacobprecon'
         [R,permR] = feval(pcmtx,H,pcoptions,DM,DG,varargin{:});
      otherwise
         error('Invalid string used for PRECONFLAG argument to TRDOG');
      end
      % We now pass kmax in from calling function
      %kmax = max(1,floor(n/2));
      if tol <= 0,
         tol = .1;
      end

      [v1,posdef,pcgit] = pcgr(DM,DG,grad,kmax,tol,...
         mtxmpy,H,R,permR,preconflag,varargin{:});
   end
   if norm(v1) > 0
      v1 = v1/norm(v1);
   end
   Z(:,1) = v1;
   if n > 1
      if (posdef < 1)
         v2 = D*sign(grad);
         if norm(v2) > 0
            v2 = v2/norm(v2);
         end
         v2 = v2 - v1*(v1'*v2);
         nrmv2 = norm(v2);
         if nrmv2 > tol2
            v2 = v2/nrmv2;
            Z(:,2) = v2;
         end
      else
         if norm(grad) > 0
            v2 = grad/norm(grad);
         else
            v2 = grad;
         end
         v2 = v2 - v1*(v1'*v2);
         nrmv2 = norm(v2);
         if nrmv2 > tol2
            v2 = v2/nrmv2;
            Z(:,2) = v2;
         end
      end
   end
end

%  REDUCE TO THE CHOSEN SUBSPACE
W = DM*Z;
switch preconflag
case 'hessprecon'
   WW = feval(mtxmpy,H,W,varargin{:});
case 'jacobprecon'
   WW = feval(mtxmpy,H,W,0,varargin{:});
otherwise
   error('Invalid string used for PRECONFLAG argument to TRDOG');
end

W = DM*WW;
MM = full(Z'*W + Z'*DG*Z);
rhs=full(Z'*grad);

%  Determine 2-D TR soln
[st,qpval,po,fcnt,lambda] = trust(rhs,MM,delta);
ss = Z*st;
s = abs(diag(D)).*ss;
s = full(s);
ssave = s;
sssave = ss;
stsave = st;

% Truncate the TR solution?
arg = (abs(s) > 0);
if isnan(s)
   error('Trust region step contains NaN''s.')
end
% No truncation if s is zero length
if isempty(find(arg))
   alpha = 1;
   mmdis = 1;
else
   mdis = inf;
   dis = max((u(arg)-x(arg))./s(arg), (l(arg)-x(arg))./s(arg));
   [mmdis,ipt] = min(dis);
   mdis = theta*mmdis;
   alpha = min(1,mdis);
end
s = alpha*s;
st = alpha*st;
ss = full(alpha*ss);
qpval1 = rhs'*st + (.5*st)'*MM*st;
if n > 1
   %   Evaluate along the reflected direction?
   qpval3 = inf;
   ssssave = mmdis*sssave;
   if norm(ssssave) < .9*delta
      r = mmdis*ssave;
      ns = ssave;
      ns(ipt) = -ns(ipt);
      nx = x+r;
      stsave = mmdis*stsave;
      qpval0 = rhs'*stsave + (.5*stsave)'*MM*stsave;
      switch preconflag
      case 'hessprecon'
         ng = feval(mtxmpy,H,r,varargin{:});
      case 'jacobprecon'
         ng = feval(mtxmpy,H,r,0,varargin{:});
      otherwise
         error('Invalid string used for PRECONFLAG argument to TRDOG');
      end

      ng = ng + g;
      ngrad = D*ng;
      ngrad = ngrad + DG*ssssave;

      %      nss is the reflected direction
      nss = sssave;
      nss(ipt) = -nss(ipt);
      ZZ(:,1) = nss/norm(nss);
      W = DM*ZZ;

      switch preconflag
      case 'hessprecon'
         WW = feval(mtxmpy,H,W,varargin{:});
      case 'jacobprecon'
         WW = feval(mtxmpy,H,W,0,varargin{:});
      otherwise
         error('Invalid string used for PRECONFLAG argument to TRDOG');
      end


      W = DM*WW;
      MM = full(ZZ'*W + ZZ'*DG*ZZ);
      nrhs=full(ZZ'*ngrad);
      [nss,tau] = quad1d(nss,ssssave,delta);
      nst = tau/norm(nss);
      ns = abs(diag(D)).*nss;
      ns = full(ns);

      %      Truncate the reflected direction?
      arg = (abs(ns) > 0);
      if isnan(ns)
         error('Reflected trust region step contains NaN''s.')
      end
      % No truncation if s is zero length
      if isempty(find(arg))
         alpha = 1;
      else
         mdis = inf;
         dis = max((u(arg)-nx(arg))./ns(arg), (l(arg)-nx(arg))./ns(arg));
         mdis = min(dis);
         mdis = theta*mdis;
         alpha = min(1,mdis);
      end
      ns = alpha*ns;
      nst = alpha*nst;
      nss = full(alpha*nss);
      qpval3 = qpval0 +  nrhs'*nst + (.5*nst)'*MM*nst;
   end

   %   Evaluate along gradient direction
   ZZ(:,1) = grad/norm(grad);
   W = DM*ZZ;

   switch preconflag
   case 'hessprecon'
      WW = feval(mtxmpy,H,W,varargin{:});
   case 'jacobprecon'
      WW = feval(mtxmpy,H,W,0,varargin{:});
   otherwise
      error('Invalid string used for PRECONFLAG argument to TRDOG');
   end


   W = DM*WW;
   MM = full(ZZ'*W + ZZ'*DG*ZZ);
   rhs=full(ZZ'*grad);
   [st,qpval,po,fcnt,lambda] = trust(rhs,MM,delta);
   ssg = ZZ*st;
   sg = abs(diag(D)).*ssg;
   sg = full(sg);

   %   Truncate the gradient direction?
   arg = (abs(sg) > 0);
   if isnan(sg)
      % No truncation if s is zero length
      error('Gradient step contains NaN''s.')
   end
   if  isempty(find(arg))
      alpha = 1;
   else
      mdis = inf;
      dis = max((u(arg)-x(arg))./sg(arg), (l(arg)-x(arg))./sg(arg));
      mdis = min(dis);
      mdis = theta*mdis;
      alpha = min(1,mdis);
   end
   sg = alpha*sg;
   st = alpha*st;
   ssg = full(alpha*ssg);
   qpval2 = rhs'*st + (.5*st)'*MM*st;
end

% Choose the best of s, sg, ns.
if qpval2 <= min(qpval1,qpval3)
   qpval = qpval2;
   s = sg;
   snod = ssg;
elseif qpval1 <= min(qpval2,qpval3)
   qpval = qpval1;
   snod = ss;
else
   qpval = qpval3;
   s = ns + r;
   snod = nss + ssssave;
end

%-----------------------------------------------------------
function[nx,tau] = quad1d(x,ss,delta)
%QUAD1D	1D quadratic zero finder for trust region step
%
% [nx,tau] = quad1d(x,ss,delta) tau is min(1,step-to-zero)
% of a 1-D quadratic ay^2 + b*y + c.
% a = x'*x; b = 2*(ss'*x); c = ss'*ss-delta^2). nx is the
% new x value, nx = tau*x;

% Algorithm:
% numer = -(b + sign(b)*sqrt(b^2-4*a*c));
% root1 = numer/(2*a);
% root2 = c/(a*root1);   % because root2*root1 = (c/a);

a = x'*x;
b = 2*(ss'*x);
c = ss'*ss-delta^2;

numer = -(b + sign(b)*sqrt(b^2-4*a*c));
warnstate = warning('off'); % Avoid divide by zero warnings
r1 = numer/(2*a);
r2 = c/(a*r1);
warning(warnstate);

tau = max(r1,r2);
tau = min(1,tau);
if tau <= 0,
   error('square root error in trdog/quad1d');
end
nx = tau*x;


%===== pcgr.m ====================================================


function[p,posdef,k] = pcgr(DM,DG,g,kmax,tol,mtxmpy,H,R,pR,callerflag,varargin);
%PCGR	Preconditioned conjugate gradients
%
% [p,posdef,k] = PCGR(DM,DG,g,kmax,tol,mtxmpy,H,R,pR) apply
% a preconditioned conjugate gradient procedure to the quadratic
%
%         q(p) = .5p'Mp + g'p, where
%
% M = DM*H*DM + DG. kmax is a bound on the number of permitted
% CG-iterations, tol is a stopping tolerance on the residual (default
% is tol = .1), mtxmpy is the function that computes products
% with the Hessian matrix H,
% and R is the cholesky factor of the preconditioner (transpose) of
% M. So, R'R approximates M(pR,pR), where pR is a permutation vector.
% On ouput p is the computed direction, posdef = 1 implies
% only positive curvature (in M) has been detected; posdef = 0
% implies p is a direction of negative curvature (for M).
% Output parameter k is the number of CG-iterations used (which
% corresponds to the number of multiplications with H).
%

% Initializations.
n = length(DG);
r = -g;
p = zeros(n,1);
val = 0;
m = 0;

% Precondition .
z = preproj(r,R,pR);
znrm = norm(z);
stoptol = tol*znrm;
inner2 = 0;
inner1 = r'*z;
posdef = 1;

kmax = max(kmax,1);  % kmax must be at least 1
% PRIMARY LOOP.
for k = 1:kmax
   if k==1
      d = z;
   else
      beta = inner1/inner2;
      d = z + beta*d;
   end
   ww = DM*d;
   switch callerflag
   case 'hessprecon'
      w = feval(mtxmpy,H,ww,varargin{:});
   case 'jacobprecon'
      w = feval(mtxmpy,H,ww,0,varargin{:});
   otherwise
      error('PCGR does not recognize this calling function.')
   end
   ww = DM*w +DG*d;
   denom = d'*ww;
   if denom <= 0
      if norm(d) == 0
         p = d;
      else
         p = d/norm(d);
      end
      posdef = 0;
      break
   else
      alpha = inner1/denom;
      p = p + alpha*d;
      r = r - alpha*ww;
   end
   z = preproj(r,R,pR);

   % Exit?
   if norm(z) <= stoptol
      break;
   end
   inner2 = inner1;
   inner1 = r'*z;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = preproj(r,RPCMTX,ppvec);
%PREPROJ Apply preconditioner
%
% w = preproj(r,RPCMTX,ppvec) Apply a preconditioner to vector r.
% The conceptual preconditioner is H(ppvec,ppvec) = RPCMTX'*RPCMTX

% Initialization
n = length(r);
one = ones(n,1);
vtol = 100*eps*one;

if nargin < 3 | isempty(ppvec)
   ppvec = (1:n);
   if nargin < 2 | isempty(RPCMTX)
      RPCMTX = speye(n);
   end
end

%  Precondition
wbar = RPCMTX'\r(ppvec);
w(ppvec,1) = RPCMTX\wbar;


%===== hmult.m ===================================================

function W = hmult(Hinfo,Y,varargin);
%HMULT	Hessian-matrix product
%
% W = HMULT(Y,Hinfo) An example of a Hessian-matrix product function
% file, e.g. Hinfo is the actual Hessian and so W = Hinfo*Y.
%
% Note: varargin is not used but must be provided in case
% the objective function has additional problem dependent
% parameters (which will be passed to this routine as well).

W = Hinfo*Y;


%===== trust.m ===================================================

function [s,val,posdef,count,lambda] = trust(g,H,delta)
%TRUST	Exact soln of trust region problem
%
% [s,val,posdef,count,lambda] = TRUST(g,H,delta) Solves the trust region
% problem: min{g^Ts + 1/2 s^THs: ||s|| <= delta}. The full
% eigen-decomposition is used; based on the secular equation,
% 1/delta - 1/(||s||) = 0. The solution is s, the value
% of the quadratic at the solution is val; posdef = 1 if H
% is pos. definite; otherwise posdef = 0. The number of
% evaluations of the secular equation is count, lambda
% is the value of the corresponding Lagrange multiplier.
%
%
% TRUST is meant to be applied to very small dimensional problems.

% INITIALIZATION
tol = 10^(-12);
tol2 = 10^(-8);
key = 0;
itbnd = 50;
lambda = 0;
n = length(g);
coeff(1:n,1) = zeros(n,1);
H = full(H);
[V,D] = eig(H);
count = 0;
eigval = diag(D);
[mineig,jmin] = min(eigval);
alpha = -V'*g;
sig = sign(alpha(jmin)) + (alpha(jmin)==0);

% POSITIVE DEFINITE CASE
if mineig > 0
   coeff = alpha ./ eigval;
   lambda = 0;
   s = V*coeff;
   posdef = 1;
   nrms = norm(s);
   if nrms <= 1.2*delta
      key = 1;
   else
      laminit = 0;
   end
else
   laminit = -mineig;
   posdef = 0;
end

% INDEFINITE CASE
if key == 0
   if seceqn(laminit,eigval,alpha,delta) > 0
      [b,c,count] = rfzero('seceqn',laminit,itbnd,eigval,alpha,delta,tol);
      vval = abs(seceqn(b,eigval,alpha,delta));
      if abs(seceqn(b,eigval,alpha,delta)) <= tol2
         lambda = b;
         key = 2;
         lam = lambda*ones(n,1);
         w = eigval + lam;
         arg1 = (w==0) & (alpha == 0);
         arg2 = (w==0) & (alpha ~= 0);
         coeff(w ~=0) = alpha(w ~=0) ./ w(w~=0);
         coeff(arg1) = 0;
         coeff(arg2) = Inf;
         coeff(isnan(coeff))=0;
         s = V*coeff;
         nrms = norm(s);
         if (nrms > 1.2*delta) | (nrms < .8*delta)
            key = 5;
            lambda = -mineig;
         end
      else
         lambda = -mineig;
         key = 3;
      end
   else
      lambda = -mineig;
      key = 4;
   end
   lam = lambda*ones(n,1);
   if (key > 2)
      arg = abs(eigval + lam) < 10 * eps * max(abs(eigval),ones(n,1));
      alpha(arg) = 0;
   end
   w = eigval + lam;
   arg1 = (w==0) & (alpha == 0); arg2 = (w==0) & (alpha ~= 0);
   coeff(w~=0) = alpha(w~=0) ./ w(w~=0);
   coeff(arg1) = zeros(length(arg1(arg1>0)),1);
   coeff(arg2) = Inf *ones(length(arg2(arg2>0)),1);
   coeff(coeff==NaN)=zeros(length(coeff(coeff==NaN)),1);
   coeff(coeff==NaN)=zeros(length(coeff(coeff==NaN)),1);
   s = V*coeff; nrms = norm(s);
   if (key > 2) & (nrms < .8*delta)
      beta = sqrt(delta^2 - nrms^2);
      s = s + beta*sig*V(:,jmin);
   end
   if (key > 2) & (nrms > 1.2*delta)
      [b,c,count] = rfzero('seceqn',laminit,itbnd,eigval,alpha,delta,tol);
      lambda = b; lam = lambda*(ones(n,1));
      w = eigval + lam;
      arg1 = (w==0) & (alpha == 0); arg2 = (w==0) & (alpha ~= 0);
      coeff(w~=0) = alpha(w~=0) ./ w(w~=0);
      coeff(arg1) = zeros(length(arg1(arg1>0)),1);
      coeff(arg2) = Inf *ones(length(arg2(arg2>0)),1);
      coeff(coeff==NaN)=zeros(length(coeff(coeff==NaN)),1);
      s = V*coeff; nrms = norm(s);
   end
end
val = g'*s + (.5*s)'*(H*s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[value] = seceqn(lambda,eigval,alpha,delta);
%SEC	Secular equation
%
% value = SEC(lambda,eigval,alpha,delta) returns the value
% of the secular equation at a set of m points lambda
%
%


%
m = length(lambda); n = length(eigval);
unn = ones(n,1); unm = ones(m,1);
M = eigval*unm' + unn*lambda'; MC = M;
MM = alpha*unm';
M(M~=0) = MM(M~=0) ./ M(M~=0);
M(MC==0) = Inf*ones(size(MC(MC==0)));
M = M.*M;
value = sqrt(unm ./ (M'*unn));
value(value==NaN) = zeros(length(value(value==NaN)),1);
value = (1/delta)*unm - value;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [b,c,itfun] = rfzero(FunFcn,x,itbnd,eigval,alpha,delta,tol,trace)
%RFZERO Find zero to the right
%
%	[b,c,itfun] = rfzero(FunFcn,x,itbnd,eigval,alpha,delta,tol,trace)
%	Zero of a function of one variable to the RIGHT of the
%       starting point x. A small modification of the M-file fzero,
%       described below, to ensure a zero to the Right of x is
%       searched for.
%
%	RFZERO is a slightly modified version of function FZERO




%	FZERO(F,X) finds a zero of f(x).  F is a string containing the
%	name of a real-valued function of a single real variable.   X is
%	a starting guess.  The value returned is near a point where F
%	changes sign.  For example, FZERO('sin',3) is pi.  Note the
%	quotes around sin.  Ordinarily, functions are defined in M-files.
%
%	An optional third argument sets the relative tolerance for the
%	convergence test.   The presence of an nonzero optional fourth
%	argument triggers a printing trace of the steps.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	C.B. Moler 1-19-86
%	Revised CBM 3-25-87, LS 12-01-88.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This algorithm was originated by T. Dekker.  An Algol 60 version,
%  with some improvements, is given by Richard Brent in "Algorithms for
%  Minimization Without Derivatives", Prentice-Hall, 1973.  A Fortran
%  version is in Forsythe, Malcolm and Moler, "Computer Methods
%  for Mathematical Computations", Prentice-Hall, 1976.
%
% Initialization
if nargin < 7, trace = 0; tol = eps; end
if nargin == 7, trace = 0; end
if trace, clc, end
itfun = 0;
%
%if x ~= 0, dx = x/20;
%if x ~= 0, dx = abs(x)/20;
if x~= 0, dx = abs(x)/2;
   %
   %else, dx = 1/20;
else, dx = 1/2;
end
%
%a = x - dx;  fa = feval(FunFcn,a,eigval,alpha,delta);
a = x; c = a;  fa = feval(FunFcn,a,eigval,alpha,delta);
itfun = itfun+1;
%
if trace, home, init = [a fa], end
b = x + dx;
b = x + 1;
fb = feval(FunFcn,b,eigval,alpha,delta);
itfun = itfun+1;
if trace, home, init = [b fb], end

% Find change of sign.

while (fa > 0) == (fb > 0)
   dx = 2*dx;
   %
   %  a = x - dx;  fa = feval(FunFcn,a);
   %  if trace, home, sign = [a fa], end
   %
   if (fa > 0) ~= (fb > 0), break, end
   b = x + dx;  fb = feval(FunFcn,b,eigval,alpha,delta);
   itfun = itfun+1;
   if trace, home, sign = [b fb], end
   if itfun > itbnd, break; end
end

fc = fb;
% Main loop, exit from middle of the loop
while fb ~= 0
   % Insure that b is the best result so far, a is the previous
   % value of b, and c is on the opposite of the zero from b.
   if (fb > 0) == (fc > 0)
      c = a;  fc = fa;
      d = b - a;  e = d;
   end
   if abs(fc) < abs(fb)
      a = b;    b = c;    c = a;
      fa = fb;  fb = fc;  fc = fa;
   end

   % Convergence test and possible exit
   %
   if itfun > itbnd, break; end
   m = 0.5*(c - b);
   toler = 2.0*tol*max(abs(b),1.0);
   if (abs(m) <= toler) + (fb == 0.0), break, end

   % Choose bisection or interpolation
   if (abs(e) < toler) + (abs(fa) <= abs(fb))
      % Bisection
      d = m;  e = m;
   else
      % Interpolation
      s = fb/fa;
      if (a == c)
         % Linear interpolation
         p = 2.0*m*s;
         q = 1.0 - s;
      else
         % Inverse quadratic interpolation
         q = fa/fc;
         r = fb/fc;
         p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
         q = (q - 1.0)*(r - 1.0)*(s - 1.0);
      end;
      if p > 0, q = -q; else p = -p; end;
      % Is interpolated point acceptable
      if (2.0*p < 3.0*m*q - abs(toler*q)) & (p < abs(0.5*e*q))
         e = d;  d = p/q;
      else
         d = m;  e = m;
      end;
   end % Interpolation

   % Next point
   a = b;
   fa = fb;
   if abs(d) > toler, b = b + d;
   else if b > c, b = b - toler;
      else b = b + toler;
      end
   end
   fb = feval(FunFcn,b,eigval,alpha,delta);
   itfun = itfun + 1;
   if trace, home, step = [b fb], end
end % Main loop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%===== perturb.m =================================================

function[pert,x,y] = perturb(x,l,u,del,y,sigma)
%PERTURB Perturb point from bounds
%
%   [PERT,X] = PERTURB(X,L,U,DEL) perturbs the current point X
%   slightly to shake it loose from tight (less than DEL away)
%   bounds U and L to be strictly feasible.
%   Called by SNLS and SFMINBX.
%
%   [PERT,X,Y] = PERTURB(X,L,U,DEL,Y,SIGMA) also perturbs the
%   reflected point Y with respect to SIGMA,

if nargin < 4
   del = 100*eps;
end

if (min(abs(u-x)) < del) | (min(abs(x-l)) < del)
   upperi = (u-x) < del;
   loweri = (x-l) < del;
   x(upperi) = x(upperi) - del;
   x(loweri) = x(loweri) + del;
   if nargin > 4
      y(upperi) = y(upperi) - del*sigma(upperi);
      y(loweri) = y(loweri) + del*sigma(loweri);
   end
   pert = 1;
else
   pert = 0;
end


%===== startx.m =================================================

function xstart = startx(u,l);
%STARTX	Box-centered point
%
% xstart = STARTX(u,l) returns centered point.

n = length(u);
onen = ones(n,1);
arg = (u > 1e12);
u(arg) = inf*onen(arg);
xstart = zeros(n,1);
arg1 = (u<inf)&(l==-inf); arg2 = (u== inf)&(l > -inf);
arg3 = (u<inf)&(l>-inf);  arg4 = (u==inf)&(l==-inf);
%
w = max(abs(u),ones(n,1));
xstart(arg1) = u(arg1) - .5*w(arg1);
%
ww = max(abs(l),ones(n,1));
xstart(arg2) = l(arg2) + .5*ww(arg2);
%
xstart(arg3)=(u(arg3)+l(arg3))/2;
xstart(arg4)=ones(length(arg4(arg4>0)),1);


%===== color.m ==================================================

function [group] = color(J,p);
%COLOR	Column partition for sparse finite differences.
%
%	 GROUP = COLOR(J,P) returns a partition of the
%   column corresponding to a coloring of the column-
%   intersection graph. GROUP(I) = J means column I is
%   colored J.
%
%	 All columns belonging to a color can be estimated
%   in a single finite difference.
%

%
[m,n] = size(J);
if nargin < 2,
   p = 1:n;
end
J = J(:,p);
group = zeros(n,1);
ncol = 0;
J = spones(J);
while any(group==0)
   % Build group for ncol
   ncol = ncol + 1;
   rows = zeros(m,1);
   index = find(group == 0);
   lenindex = length(index);
   for i = 1:lenindex
      k = index(i);
      inner = J(:,k)'*rows;
      if inner == 0
         group(k) = ncol;
         rows = rows + J(:,k);
      end
   end
end
group(p)= group;


%===== sfd.m ====================================================

function[H] = sfd(x,grad,H,group,alpha,funfcn,varargin)
%SFD	Sparse Hessian via finite gradient differences
%
% H = sfd(x,grad,H,group,fdata,fun) returns the
% sparse finite difference approximation H of a Hessian matrix
% of function 'fun'  at current point x.
% Vector group indicates how to use sparse finite differencing:
% group(i) = j means that column i belongs to group (or color) j.
% Each group (or color) corresponds to a finite gradient difference.
% fdata is a data array (possibly) needed by function 'fun'.
%
% H = sfd(x,grad,H,group,fdata,fun,alpha) overrides the default
% finite differencing stepsize.
%

xcurr = x(:);  % Preserve x so we know what funfcn expects
scalealpha = 0;
[m,n] = size(H);
v = zeros(n,1);
ncol = max(group); epsi = sqrt(eps);
if isempty(alpha)
    scalealpha = 1;
    alpha = ones(ncol,1)*sqrt(eps);
end
H = spones(H); d = zeros(n,1);
for k = 1:ncol
    d = (group == k);
    if scalealpha
        xnrm = norm(xcurr(d));
        xnrm = max(xnrm,1);
        alpha(k) = alpha(k)*xnrm;
    end
    y = xcurr + alpha(k)*d;

    % Make x conform to user-x
    x(:) = y;
    %[dummy,v] = feval(fun,y);
    switch funfcn{1}
    case 'fun'
        error('should not reach this')
    case 'fungrad'
        [dummy,v(:)] = feval(funfcn{3},x,varargin{:});
        %OPTIONS(11)=OPTIONS(11)+1;
    case 'fun_then_grad'
        %newval = feval(funfcn{3},x,varargin{:});
        v(:) = feval(funfcn{4},x,varargin{:});
        % OPTIONS(11)=OPTIONS(11)+1;
    otherwise
        error('Undefined calltype in FMINUNC');
    end


    w = (v-grad)/alpha(k);
    cols = find(d);
    lpoint = length(cols);
    A = sparse(m,n);
    A(:,cols) = H(:,cols);
    H(:,cols) = H(:,cols) - A(:,cols);
    [i,j,val] = find(A);
    [p,ind] = sort(i);
    val(ind) = w(p);
    A = sparse(i,j,full(val),m,n);
    H = H + A;
end
H = (H+H')/2;   % symmetricize


%===== hprecon.m ================================================

function[R,pvec] = hprecon(H,upperbandw,DM,DG,varargin);
%HPRECON Sparse Cholesky factor of H-preconditioner
%
%   [R,PVEC] = HPRECON(H,UPPERBANDW,DM,DG) computes the
%   sparse Cholesky factor (transpose of a (usually) banded
%   preconditioner of square matrix M
%                M = DM*H*DM + DG
%   wehre DM and DG are non-negative sparse diagonal matrices.
%   R'*R approximates M(pvec,pvec), i.e.
%          R'*R = M(pvec,pvec)
%
%   H may not be the true Hessian.  If H is the same size as the
%   true Hessian, H will be used in computing the preconditioner R.
%   Otherwise, compute a diagonal preconditioner for
%               M = DM*DM + DG
%
%   If 0 < UPPERBANDW <  n then the upper bandwidth of
%   R is UPPERBANDW. If UPPERBANDW >= n then the structure of R
%   corresponds to a sparse Cholesky factorization of H
%   using the symmmd ordering (the ordering is returned in PVEC).
%

%   Default preconditioner for SFMINBX and SQPMIN.

if nargin < 1,
   error('hprecon requires at least 1 input parameter.');
end
if nargin <2,
    upperbandw = 0;
    if nargin < 3
        DM = [];
        if nargin < 4
            DG = [];
        end, end, end

[rows,cols] = size(H);
n = length(DM);
% In case "H" isn't really H, but something else to use with HessMult function.
if ~isnumeric(H) | ~isequal(n,rows) | ~isequal(n,cols)
    % H is not the right size; ignore requested bandwidth and compute
    % diagonal preconditioner based only on DM and DG.
    pvec = (1:n);
    d1 = full(diag(DM));  % full vector
    d2 = full(diag(DG));
    dd = sqrt(d1.*d1 + abs(d2));
    R = sparse(1:n,1:n,dd);
    return
end

H = DM*H*DM + DG;
pvec = (1:n);
epsi = .0001*ones(n,1);
info = 1;

if upperbandw >= n-1 % Try complete approximation to H
   pvec = symmmd(H);
   ddiag = diag(H);
   mind = min(ddiag);
   lambda = 0;
   if mind < 0,
      lambda = -mind + .001;
   end
   while info > 0
      H = H + lambda*speye(n);
      [R,info] = chol(H(pvec,pvec));
      lambda = lambda + 10;
   end
elseif (upperbandw > 0) & ( upperbandw < n-1) % Banded approximation to H
   % Banded approximation
   lambda = 0;
   ddiag = diag(H);
   mind = min(ddiag);
   if mind < 0,
      lambda = -mind + .001;
   end
   H = tril(triu(H,-upperbandw),upperbandw);
   while info > 0
      H = H + lambda*speye(n);
      [R,info] = chol(H);
      lambda = 4*lambda;
      if lambda <= .001,
         lambda = 1;
      end
   end
elseif upperbandw == 0 % diagonal approximation for H
   dnrms = sqrt(sum(H.*H))';
   d = max(sqrt(dnrms),epsi);
   R = sparse(1:n,1:n,full(d));
   pvec = (1:n);
else
    error('upperbandw must be >= 0.')
end
