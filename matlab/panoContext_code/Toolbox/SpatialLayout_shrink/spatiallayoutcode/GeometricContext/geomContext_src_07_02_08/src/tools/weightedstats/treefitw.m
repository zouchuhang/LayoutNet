function Tree=treefitw(X,y,w, equivsample, varargin)
%TREEFIT Fit a tree-based model for classification or regression.
%   T = TREEFIT(X,Y) creates a decision tree T for predicting response Y
%   as a function of predictors X.  X is an N-by-M matrix of predictor
%   values.  Y is either a vector of N response values (for regression),
%   or a character array or cell array of strings containing N class
%   names (for classification).  Either way, T is binary tree where each
%   non-terminal node is split based on the values of a column of X.  NaN
%   values in X or Y are taken to be missing values, and observations with
%   any missing values are not used in the fit.
%
%   W is the weight vector with sum(W) = N
%   Equivsample is the equivalent sample size (add that number to each node
%   for each class.
%
%   T = TREEFIT(X,Y,W,0,'PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%
%   For all trees:
%      'catidx'     Vector of indices of the columns of X that are to be
%                   treated as unordered categorical variables
%      'method'     Either 'classification' (default if Y is text) or
%                   'regression' (default if Y is numeric)
%      'splitmin'   A number N such that impure nodes must have N or more
%                   observations to be split (default 10)
%      'prune'      'on' (default) to compute the full tree and the optimal
%                   sequence of pruned subtrees, or 'off' for the full tree
%                   without pruning
%      'maxnodes'   The maximum number of splitting nodes in the tree
%
%   For classification trees only:
%      'cost'       Square matrix C, C(i,j) is the cost of classifying
%                   a point into class j if its true class is i (default
%                   has C(i,j)=1 if i~=j, and C(i,j)=0 if i=j).  Alternatively
%                   this value can be a structure S having two fields:  S.group
%                   continaing the group names as a character array or cell
%                   array of strings, and S.cost containing the cost matrix C.
%      'splitcriterion'  Criterion for choosing a split, either 'gdi' (default)
%                   for Gini's diversity index, 'twoing' for the twoing rule,
%                   or 'deviance' for maximum deviance reduction
%      'priorprob'  Prior probabilities for each class, specified as a
%                   vector (one value for each distinct group name) or as a
%                   structure S with two fields:  S.group containing the group
%                   names as a character array or cell array of strings, and
%                   S.prob containing a a vector of corresponding probabilities
%
%   Example:  Create classification tree for Fisher's iris data.
%      load fisheriris;
%      t = treefit(meas, species);
%      treedisp(t,'names',{'SL' 'SW' 'PL' 'PW'});
%
%   See also TREEDISP, TREEPRUNE, TREETEST, TREEVAL.

%   Reference:  Breiman et al. (1993), "Classification and Regression
%   Trees," Chapman and Hall, Boca Raton.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 1.1 $  $Date: 2004/08/23 14:32:11 $

% Process inputs
if isnumeric(y)
   Method = 'regression';
else
   Method = 'classification';
   if (ischar(y))
      y = cellstr(y);
   end
end
if ~isnumeric(X)
   error('X must be a numeric matrix.');
end
okargs =   {'priorprob' 'cost' 'splitcriterion' 'splitmin' 'catidx' 'prune' 'method' 'maxnodes'};
defaults = {[]          []     'gdi'            10         []       'on'    Method  1000};
[emsg Prior Cost Criterion Splitmin Catidx Prune Method maxNodes] = ...
                statgetargs(okargs,defaults,varargin{:});
error(emsg);

if ~isstr(Method) | isempty(Method) | ~(Method(1)=='c' | Method(1)=='r')
   error('Value of ''method'' parameter must be ''classification'' or ''regression''.');
elseif Method(1)=='c'
   Method = 'classification';
else
   Method = 'regression';
end

t = any(isnan(X),2);
if isequal(Method,'regression')
   t = t | isnan(y);
end
if any(t)
   disp(['nan: ' num2str(find(t==1)')])
   X(t,:) = [];
   y(t) = [];
end

[N,nvars] = size(X);
doclass = isequal(Method(1),'c');
if doclass
   switch(Criterion)
    %                Criterion function   Is it an impurity measure?
    %                ------------------   --------------------------
    case 'gdi',      critfun = @gdi;      impurity = 1;
    case 'twoing',   critfun = @twoing;   impurity = 0;
    case 'deviance', critfun = @deviance; impurity = 1;
    otherwise,     error('Bad value for ''splitcriterion'' parameter.')
   end
   
   % Get binary matrix, C(i,j)==1 means point i is in class j
   if islogical(y)
      y = double(y);
   end
   [y,cnames] = grp2idx(y);   % find groups only after NaNs removed from X
   if any(isnan(y))
      t = isnan(y);
      y(t) = [];
      X(t,:) = [];
      N = size(X,1);
   end
   nclasses = max(y);
   C = zeros(N,nclasses);
   C(sub2ind([N nclasses],(1:N)',y)) = 1;   
   for cindex = 1:nclasses
       C(:, cindex) =C(:, cindex).*w;
   end
   Nj = sum(C,1);
else
   C = y(:);
end

% Tree structure fields ([C] only for classification trees):
%  .method     method
%  .node       node number
%  .parent     parent node number
%  .class      class assignment for points in this node if treated as a leaf
%  .var        column j of X matrix to be split, or 0 for a leaf node,
%              or -j to treat column j as categorical
%  .cut        cutoff value for split (Xj<cutoff goes to left child node),
%              or index into catsplit if var is negative
%  .children   matrix of child nodes (2 cols, 1st is left child)
%  .nodeprob   probability p(t) for this node
%  .nodeerr    resubstitution error estimate r(t) for this node
%  .risk       R(t) = p(t)*r(t)
%  .nodesize   number of points at this node
%  .prunelist  list of indices that define pruned subtrees.  One entry per
%              node.  If prunelist(j)=k then, at the kth level of pruning,
%              the jth node becomes a leaf (or drops off the tree if its
%              parent also gets pruned).
%  .alpha      vector of complexity parameters for each pruning cut
%  .ntermnodes vector of terminal node counts for each pruning cut
%  .catsplit   call array for categorical splits,
%              left categories in column 1 and right categories in column 2
%  .classprob  [C] vector of class probabilities
%  .classname  [C] names of each class
%  .classcount [C] count of members of each class
%  .nclasses   [C] number of classes

nodenumber = zeros(N,1);
parent = zeros(N,1);
yfitnode = zeros(N,1);
cutvar = zeros(N,1);
cutpoint = zeros(N,1);
children = zeros(N,2);
nodeprob = zeros(N,1);
resuberr = zeros(N,1);
risk = zeros(N,1);
nodesize = zeros(N,1);
if doclass
   classprob = zeros(N,nclasses);
   classcount = zeros(N,nclasses);
end
catsplit = cell(0,2);
iscat = zeros(nvars,1); iscat(Catidx) = 1;

nodenumber(1) = 1;

assignednode = ones(N,1);
nextunusednode = 2;

if doclass
   % Get default or specified prior class probabilities
   Prior = Prior(:)';
   haveprior = true;
   if isempty(Prior)
      Prior = Nj / N;
      haveprior = false;
   elseif isequal(Prior,'equal')
      Prior = ones(1,nclasses) / nclasses;

   elseif isstruct(Prior)
      if ~isfield(Prior,'group') | ~isfield(Prior,'prob')
         error('Missing field in structure value for ''priorprob'' parameter.');
      end
      idx = getclassindex(cnames,Prior.group);
      if any(idx==0)
         j = find(idx==0);
         error(sprintf('Missing prior probability for group ''%s''.',...
                       cnames{j(1)}));
      end
      Prior = Prior.prob(idx);
   end
   if length(Prior)~=nclasses | any(Prior<0) | sum(Prior)==0 ...
                              | ~isnumeric(Prior)
      error(sprintf(...
       'Value of ''priorprob'' parameter must be a vector of %d probabilities.',...
        nclasses));
   else
      Prior = Prior / sum(Prior);
   end

   % Get default or specified misclassification costs
   havecosts = true;
   if isempty(Cost)
      Cost = ones(nclasses) - eye(nclasses);
      havecosts = false;
   else
      if isstruct(Cost)
         if ~isfield(Cost,'group') | ~isfield(Cost,'cost')
            error('Missing field in structure value for ''cost'' parameter.');
         end
         idx = getclassindex(cnames,Cost.group);
         if any(idx==0)
            j = find(idx==0);
            error(sprintf('Missing misclassification cost for group ''%s''.',...
                          cnames{j(1)}));
         end
         Cost = Cost.cost(idx,idx);
      end
      if ~isequal(size(Cost),nclasses*ones(1,2))
         error(sprintf('Misclassification cost matrix must be %d-by-%d.',...
                       nclasses,nclasses));
      elseif any(diag(Cost)~=0)
         error('Misclassification cost matrix must have zeros on the diagonal.');
      elseif any(Cost<0)
         error('Misclassification cost matrix must contain non-negative values.');
      end
   end
   
   % Adjust priors if required to take misclassification costs into account
   adjprior = Prior;
   if havecosts
      Cj = sum(Cost,2)';
      pc = Cj .* Prior;
      adjprior = pc / sum(pc);
   end
end

% Keep processing nodes until done
tnode = 1;
while(tnode < nextunusednode)
   % Record information about this node
   noderows = find(assignednode==tnode);
   Nnode = length(noderows);
   Cnode = C(noderows,:);
   NnodeC = sum(Cnode>0, 1);
   if doclass
      % Compute class probabilities and related statistics for this node
      %Njt = sum(Cnode,1);    % number in class j at node t
      Njt = sum(Cnode,1) + equivsample;    % number in class j at node t
      Pjandt = Prior .* Njt ./ (Nj+equivsample*nclasses);
      Pjgivent = Pjandt / sum(Pjandt);
      misclasscost = Pjgivent * Cost;
      [mincost,nodeclass] = min(misclasscost);
      yfitnode(tnode) = nodeclass;
      Pt = sum(Pjandt);
      nodeprob(tnode) = Pt;
      classprob(tnode,:) = Pjgivent;
      classcount(tnode,:) = Njt;
      pratio = adjprior ./ Nj;
      % was ... impure = sum(Pjgivent>0)>1;
      impure = sum(NnodeC>2)>1; % must be "significantly" impure
   else
      % Compute variance and related statistics for this node
      ybar = mean(Cnode);
      yfitnode(tnode) = ybar;
      nodeprob(tnode) = Nnode/N;
      sst = norm(Cnode-ybar)^2;   % total sum of squares at this node
      mincost = sst / Nnode;
      impure = (mincost > 1e-6 * resuberr(1));
   end
   bestcrit          = -Inf;
   nodesize(tnode)   = Nnode;
   resuberr(tnode)   = mincost;
   risk(tnode)       = nodeprob(tnode) * resuberr(tnode);
   cutvar(tnode)     = 0;
   cutpoint(tnode)   = 0;
   children(tnode,:) = 0;
   
   % Consider splitting this node
   if Nnode>=Splitmin & impure  & tnode <= maxNodes % split only large impure nodes
      Xnode = X(noderows,:);
      bestvar = 0;
      bestcut = 0;

      % Find the best of all possible splits
      for jvar=1:nvars
         x = Xnode(:,jvar);            % get jth x variable
         [x,idx] = sort(x);            % sort it
         xcat = iscat(jvar);         
         if doclass
            Ccum = cumsum(Cnode(idx,:));         % cum. class counts
            [critval,cutval]=Ccritval(x,Ccum,xcat,pratio,Pt,impurity,critfun);
         else
            ycum = cumsum(Cnode(idx,:) - ybar);  % centered response cum. sum
            [critval,cutval]=Rcritval(x,ycum,xcat);
         end

         % Change best split if this one is best so far
         if critval>bestcrit
            bestcrit = critval;
            bestvar = jvar;
            bestcut = cutval;
         end
      end

      % Split this node using the best rule found
      if bestvar~=0
         x = Xnode(:,bestvar);
         if ~iscat(bestvar)
            cutvar(tnode) = bestvar;
            cutpoint(tnode) = bestcut;
            leftside = x<=bestcut;
            rightside = ~leftside;
         else
            cutvar(tnode) = -bestvar;          % negative indicates cat. var. split
            ncatsplit = size(catsplit,1) + 1;  % index into catsplit cell array
            cutpoint(tnode) = ncatsplit;
            catsplit(ncatsplit,:) = bestcut;
            leftside = ismember(x,bestcut{1});
            rightside = ismember(x,bestcut{2});
         end
         children(tnode,:) = nextunusednode + (0:1);
         assignednode(noderows(leftside)) = nextunusednode;
         assignednode(noderows(rightside)) = nextunusednode+1;
         nodenumber(nextunusednode+(0:1)) = nextunusednode+(0:1)';
         parent(nextunusednode+(0:1)) = tnode;
         nextunusednode = nextunusednode+2;
      end
   end
   tnode = tnode + 1;
end

topnode        = nextunusednode - 1;
Tree.method    = Method;
Tree.node      = nodenumber(1:topnode);
Tree.parent    = parent(1:topnode);
Tree.class     = yfitnode(1:topnode);
Tree.var       = cutvar(1:topnode);
Tree.cut       = cutpoint(1:topnode);
Tree.children  = children(1:topnode,:);
Tree.nodeprob  = nodeprob(1:topnode);
Tree.nodeerr   = resuberr(1:topnode);
Tree.risk      = risk(1:topnode);
Tree.nodesize  = nodesize(1:topnode);
Tree.npred     = nvars;
Tree.catcols   = Catidx;
if doclass
   if ~haveprior, Prior=[]; end
   Tree.prior     = Prior;
   Tree.nclasses  = nclasses;
   Tree.cost      = Cost;
   Tree.classprob = classprob(1:topnode,:);
   Tree.classcount= classcount(1:topnode,:);
   Tree.classname = cnames;
end

Tree.catsplit  = catsplit; % list of all categorical predictor splits

Tree = removebadsplits(Tree);

if isequal(Prune,'on')
   Tree = treeprune(Tree);
end

%----------------------------------------------------
function v=gdi(p)
%GDI Gini diversity index

v=1-sum(p.^2,2);

%----------------------------------------------------
function v=twoing(Pleft, P1, Pright, P2)
%TWOING Twoing index

v = 0.25 * Pleft .* Pright .* sum(abs(P1-P2),2).^2;

%----------------------------------------------------
function v=deviance(p)
%DEVIANCE Deviance

v = -2 * sum(p .* log(max(p,eps)), 2);

%----------------------------------------------------
function [critval,cutval]=Ccritval(x,Ccum,iscat,pratio,Pt,impurity,critfun)
%CCRITVAL Get critical value for splitting node in classification tree.
   
% First get all possible split points
Ncum = (1:length(x))';
rows = Ncum(diff(x)>0);
if isempty(rows)
   critval = -Inf;
   cutval = 0;
   return
end

% Get arrays showing left/right class membership at each split
nsplits = length(rows);
if iscat
   % A picks out all category subsets including the 1st, but not the whole set
   A = ones(2^nsplits,nsplits+1);
   A(:,2:end) = fullfact(2*ones(1,nsplits)) - 1;
   A(end,:) = [];
   
   % B contains the class counts in each category
   t = [rows; size(Ccum,1)];
   B = Ccum(t,:);
   B(2:end,:) = B(2:end,:) - B(1:end-1,:);
   
   Csplit1 = A*B;
   nsplits = size(Csplit1,1);
   allx = x(t);
else
   % Split between each pair of distinct ordered values
   Csplit1 = Ccum(rows,:);
end
Csplit2 = repmat(Ccum(end,:),nsplits,1) - Csplit1;

% Get left/right class probabilities at each split
temp = repmat(pratio,nsplits,1);
P1 = temp .* Csplit1;
P2 = temp .* Csplit2;
Ptleft  = sum(P1,2);
Ptright = sum(P2,2);
nclasses = size(P1,2);
P1 = P1 ./ max(repmat(Ptleft,1,nclasses), 1E-10); % max added by DWH
P2 = P2 ./ max(repmat(Ptright,1,nclasses), 1E-10);

% Get left/right node probabilities
Pleft = Ptleft ./ Pt;
Pright = 1 - Pleft;

% Evaluate criterion as impurity or otherwise
if impurity
   crit = - Pleft.*feval(critfun,P1) - Pright.*feval(critfun,P2);
else
   crit = feval(critfun, Pleft, P1, Pright, P2);
end

% Return best split point
critval = max(crit);
maxloc = find(crit==critval);
if length(maxloc)>1
   maxloc = maxloc(1+floor(length(maxloc)*rand));
end
if iscat
   t = logical(A(maxloc,:));
   xleft = allx(t);
   xright = allx(~t);
   cutval = {xleft' xright'};
else
   cutloc = rows(maxloc);
   cutval = (x(cutloc) + x(cutloc+1))/2;
end

%----------------------------------------------------
function [critval,cutval]=Rcritval(x,Ycum,iscat)
%RCRITVAL Get critical value for splitting node in regression tree.
   
% First get all possible split points
Ncum = (1:length(x))';
rows = Ncum(diff(x)>0);
if isempty(rows)
   critval = -Inf;
   cutval = 0;
   return
end

% Get arrays showing left/right class membership at each split
nsplits = length(rows);
if iscat
   % A picks out all category subsets including the 1st, but not the whole set
   A = ones(2^nsplits,nsplits+1);
   A(:,2:end) = fullfact(2*ones(1,nsplits)) - 1;
   A(end,:) = [];
   
   % B contains the category sums
   t = [rows; size(Ycum,1)];
   B = Ycum(t,:);
   B(2:end,:) = B(2:end,:) - B(1:end-1,:);
   
   Ysplit1 = A*B;
   n1 = A*[t(1);diff(t)];
   allx = x(t);               % take one x value from each unique set
else
   % Split between each pair of distinct ordered values
   Ysplit1 = Ycum(rows,:);
   n1 = rows;
end

% Get left/right means
N = Ncum(end);
mu1 = Ysplit1 ./ n1;
mu2 = (Ycum(end) - Ysplit1) ./ (N - n1);

ssx = n1.*mu1.^2 + (N-n1).*mu2.^2;
critval = max(ssx);
maxloc = find(ssx==critval);
if length(maxloc)>1
   maxloc = maxloc(1+floor(length(maxloc)*rand));
end
if iscat
   t = logical(A(maxloc,:));
   xleft = allx(t);
   xright = allx(~t);
   cutval = {xleft' xright'};
else
   cutloc = rows(maxloc);
   cutval = (x(cutloc) + x(cutloc+1))/2;
end

% --------------------------------------
function Tree = removebadsplits(Tree)
%REMOVEBADSPLITS Remove splits that contribute nothing to the tree.

N = length(Tree.node);
isleaf = (Tree.var==0)';   % no split variable implies leaf node
isntpruned = true(1,N);
doprune = false(1,N);
adjfactor = (1 - 100*eps);
risk = Tree.risk';

% Work up from the bottom of the tree
while(true)
   % Find "twigs" with two leaf children
   leafs = find(isleaf & isntpruned);
   branches = find(~isleaf & isntpruned);
   twig = branches(sum(isleaf(Tree.children(branches,:)),2) == 2);
   if isempty(twig)
      break;            % must have just the root node left
   end
   
   % Find twigs to "unsplit" if the error of the twig is no larger
   % than the sum of the errors of the children
   Rtwig = risk(twig);
   kids = Tree.children(twig,:);
   Rsplit = sum(risk(kids),2);
   unsplit = Rsplit >= Rtwig'*adjfactor;
   if any(unsplit)
      % Mark children as pruned, and mark twig as now a leaf
      isntpruned(kids(unsplit,:)) = 0;
      twig = twig(unsplit);   % only these to be marked on next 2 lines
      isleaf(twig) = 1;
      doprune(twig) = 1;
   else
      break;
   end
end

% Remove splits that are useless
if any(doprune)
   Tree = treeprune(Tree,'nodes',find(doprune));
end
