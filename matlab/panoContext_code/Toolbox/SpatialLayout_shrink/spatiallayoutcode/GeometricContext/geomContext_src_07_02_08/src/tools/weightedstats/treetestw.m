function [cost,secost,ntnodes,bestlevel] = treetestw(Tree,TorCorR,X,Y,w,varargin)
%TREETEST Compute error rate for tree.
%   COST = TREETEST(T,'resubstitution') computes the cost of the tree T
%   using a resubstitution method.  T is a decision tree as created by
%   the TREEFIT function.  The cost of the tree is the sum over all
%   terminal nodes of the estimated probability of that node times the
%   node's cost.  If T is a classification tree, the cost of a node is
%   the sum of the misclassification costs of the observations in
%   that node.  If T is a regression tree, the cost of a node is the
%   average squared error over the observations in that node.  COST is
%   a vector of cost values for each subtree in the optimal pruning
%   sequence for T.  The resubstitution cost is based on the same
%   sample that was used to create the original tree, so it under-
%   estimates the likely cost of applying the tree to new data.
%
%   COST = TREETEST(T,'test',X,Y) uses the predictor matrix X and
%   response Y as a test sample, applies the decision tree T to that
%   sample, and returns a vector COST of cost values computed for the
%   test sample.  X and Y should not be the same as the learning sample,
%   which is the sample that was used to fit the tree T.
%
%   COST = TREETEST(T,'crossvalidate',X,Y) uses 10-fold cross-validation to
%   compute the cost vector.  X and Y should be the learning sample, which
%   is the sample that was  used to fit the tree T.  The function
%   partitions the sample into 10 subsamples, chosen randomly but with
%   roughly equal size.  For classification trees the subsamples also have
%   roughly the same class proportions.  For each subsample, TREETEST fits
%   a tree to the remaining data and uses it to predict the subsample.  It
%   pools the information from all subsamples to compute the cost for the
%   whole sample.  w is the weight vector for the samples.  If no weights
%   the original function treetest should be used instead.  
%
%   [COST,SECOST,NTNODES,BESTLEVEL] = TREETEST(...) also returns the vector
%   SECOST containing the standard error of each COST value, the vector
%   NTNODES containing number of terminal nodes for each subtree, and the
%   scalar BESTLEVEL containing the estimated best level of pruning.
%   BESTLEVEL=0 means no pruning (i.e. the full unpruned tree).  The best
%   level is the one that produces the smallest tree that is within one
%   standard error of the minimum-cost subtree.
%
%   [...] = TREETEST(...,'PARAM1',val1,'PARAM2',val2,...) specifies
%   optional parameter name/value pairs chosen from the following:
%
%      'nsamples'   The number of cross-validation samples (default 10)
%      'treesize'   Either 'se' (the default) to choose the smallest
%                   tree whose cost is within one standard error of the
%                   minimum cost, or 'min' to choose the minimal cost tree
%                   (not meaningful for resubstitution error calculations)
%
%   Example:  Find best tree for Fisher's iris data using cross-validation.
%             The solid line shows the estimated cost for each tree size,
%             the dashed line marks 1 standard error above the minimum,
%             and the square marks the smallest tree under the dashed line.
%      % Start with a large tree
%      load fisheriris;
%      t = treefit(meas,species','splitmin',5);
%
%      % Find the minimum-cost tree
%      [c,s,n,best] = treetest(t,'cross',meas,species);
%      tmin = treeprune(t,'level',best);
%
%      % Plot smallest tree within 1 std. error of minimum cost tree
%      [mincost,minloc] = min(c);
%      plot(n,c,'b-o', n(best+1),c(best+1),'bs',...
%           n,(mincost+s(minloc))*ones(size(n)),'k--');
%      xlabel('Tree size (number of terminal nodes)')
%      ylabel('Cost')
%
%   See also TREEFIT, TREEDISP, TREEPRUNE, TREEVAL.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 1.1 $  $Date: 2004/08/23 14:32:11 $

if nargin<2, error('Not enough arguments.'); end
if ~isstruct(Tree) | ~isfield(Tree,'method')
   error('First argument must be a decision tree.');
end
if ~isstr(TorCorR) | ~(TorCorR(1)=='t' | TorCorR(1)=='c' | TorCorR=='r')
   error('Second argument must be ''test'', ''crossvalidate'', or ''resubstitution''.');
end
if TorCorR(1)=='t' & nargin<4
   error('Not enough arguments.  Need X and Y for the test sample.');
elseif TorCorR(1)=='c' & nargin<4
   error('Not enough arguments.  Need X and Y from the learning sample.');
end
doclass = isequal(Tree.method,'classification');
if TorCorR(1)~='r'
   if ~ischar(Y) & prod(size(Y))~=length(Y)
      error('Y must be a vector.');
   else
      if iscell(Y) | isnumeric(Y)
         n = length(Y);
      else
         n = size(Y,1);
      end
      if size(X,1)~=n
         error('There must be one Y value for each row of X.');
      end
   end
end

okargs =   {'nsamples' 'treesize'};
defaults = {10         'se'};
[emsg ncv treesize] = statgetargs(okargs,defaults,varargin{:});
error(emsg);

if ~isnumeric(ncv) | prod(size(ncv))~=1 | ncv<2 | ncv~=round(ncv)
   error('Value of ''nsamples'' argument must be an integer 2 or larger.');
end
if ~isstr(treesize) | ~(treesize(1)=='s' | treesize(1)=='m')
   error('Value of ''treesize'' argument must be ''se'' or ''min''.');
end

% Get complexity parameters for all pruned subtrees
if ~isfield(Tree,'alpha')
   Tree = treeprune(Tree);
end

% Remove missing values
if nargin>=4
   t = any(isnan(X),2);
   if isequal(Tree.method,'classification')
      Yold = Y;
      Y = classname2id(Y,Tree.classname);
      if any(Y==0)
         bad = find(Y==0);
         bad = Yold(bad(1));
         if isnumeric(bad)
            bad = num2str(bad);
         elseif iscell(bad)
            bad = bad{1};
         end
         error(sprintf(...
             'At least one Y value (''%s'') is incompatible with the tree.',... 
                       bad));
      end
   end
                   
   t = t | isnan(Y);
   if any(t)
      X(t,:) = [];
      Y(t,:) = [];
   end
end

% Do proper type of testing (error estimation)
switch(TorCorR(1))
 case 't', [cost,secost] = testtree(Tree,X,Y);
 case 'c', [cost,secost] = cvtree(Tree,X,Y,w,ncv);
 case 'r', [cost,secost] = resubinfo(Tree); treesize = 'm';
end

cost = cost(:);
secost = secost(:);
if nargout>=3
   ntnodes = Tree.ntermnodes(:);
end
if nargout>=4
   bestlevel = selecttree(cost,secost,treesize(1)) - 1;
end

% ---------------------------------------------------------
function [resuberr,seresub] = resubinfo(Tree)
%RESUBINFO Compute error rates for tree using resubstitution error.

% Get complexity parameters for all pruned subtrees
nsub = 1+max(Tree.prunelist);

% Get error rate for each subtree in this sequence
resuberr = zeros(nsub,1);
for j=1:nsub;
   Tj = treeprune(Tree,'level',j-1);
   leaves = Tj.node(Tj.var==0);
   resuberr(j) = sum(Tj.risk(leaves));
end
seresub = zeros(size(resuberr));

% ---------------------------------------------------------------
function [testerr,seerr] = testtree(Tree,X,id)
%TESTTREE Compute error rates for tree using test sample.
%   The id variable is the class id for classification, or the y variable
%   for regression.

% Get pruning sequence and compute fitted values for the whole sequence
nsub = 1 + max(Tree.prunelist);
yfit = treeval(Tree,X,(0:nsub-1));

doclass = isequal(Tree.method,'classification');
if doclass   % get info required for classification
   nclasses = Tree.nclasses;
   cost = Tree.cost;
   prior = Tree.prior(:);
   if isempty(prior)
      prior = Tree.classcount(1,:)' / Tree.nodesize(1);
   end
   Njtest = histc(id,1:nclasses);
   adjprior = (prior ./ max(eps,Njtest))';
end

% Compute error statistics
if doclass
   testerr = zeros(nsub,1);
   seerr = zeros(nsub,1);
   for k=nsub:-1:1;
      % M(i,j) counts class i items classified as class j
      M = full(sparse(id,yfit(:,k),1,nclasses,nclasses));
   
      % Compute loss for this classification
      loss = sum(cost .* M, 2);
      losssq = sum(cost.^2 .* M, 2);
      s2 = losssq  - loss.^2 ./ Njtest;
      
      testerr(k) = adjprior * loss;
      seerr(k) = sqrt(adjprior.^2 * s2);
   end
else
   N = size(X,1);
   E = (yfit - repmat(id,1,size(yfit,2))).^2;
   testerr = mean(E,1);
   s2 = sum((E - repmat(testerr,size(E,1),1)).^2,1) / N;
   seerr = sqrt(s2/N);
end

% ---------------------------------------------------------------
function [cverr,secverr] = cvtree(Tree,X,id,w,ncv)
%CVTREE Compute error rates for tree using cross-validaiton.
%   [CVERR,SECVERR] = CVTREE(TREE,X,ID,NCV)

% Get geometric means of the alpha boundary points
alpha = Tree.alpha;
avgalpha = [sqrt(alpha(1:end-1) .* alpha(2:end)); Inf];

% Loop over cross-validation samples
N = size(X,1);
ntrees = length(avgalpha);
cverr = zeros(ntrees,1);
secverr = zeros(ntrees,1);
cvid = 1 + mod((1:N),ncv);

doclass = isequal(Tree.method,'classification');
if doclass
   % Use a random permutation with fixed category proportions
   idrand = id + rand(size(id));
   [stdid,idx] = sort(idrand);
   cvid = cvid(idx);
   args = {'prior',Tree.prior, 'cost',Tree.cost, 'prune','on'};
else   
   % Use a random permutation with fixed numbers per cross-validation sample
   cvid = cvid(randperm(N));
   args = {'prune','on'};
end

% Get predicted values using cross-validation samples
cvclass = zeros(N,ntrees);
for j=1:ncv
   % Use jth group as a test, train on the others
   testrows = (cvid == j);
   trainrows = ~testrows;
   testsize = sum(testrows);

   % Get a sequence of pruned trees for the training set
   Tj = treefitw(X(trainrows,:),id(trainrows),w(trainrows),0,...
                'method',Tree.method, 'catidx',Tree.catcols, args{:});

   % Get classifications based on each subtree that we require
   treesneeded = findsubtree(Tj,avgalpha);
   cvclass(testrows,:) = treeval(Tj,X(testrows,:),treesneeded-1);
end

% Compute output statistics based on those predictions
if doclass
   Nj = Tree.classcount(1,:)';
   prior = Tree.prior;
   if isempty(prior)
      prior = Nj' / N;
   end
   adjprior = (prior ./ Nj');
   nclasses = length(prior);
   cost = Tree.cost;
   sz = size(cost);
   w = w / sum(w);
   for k=1:ntrees      
      loss = sum((cvclass(:, k)~=id).*w);
      losssq = sum(((cvclass(:, k)~=id).^2).*w);
      s2 = losssq - loss.^2;
      cverr(k) = loss;
      secverr(k) = sqrt(s2);
      %M = full(sparse(id,cvclass(:,k),1,nclasses,nclasses));
      %loss = sum(cost .* M, 2);
      %losssq = sum(cost.^2 .* M, 2);
      %s2 = losssq - loss.^2 ./ Nj;
      %cverr(k) = adjprior * loss;
      %secverr(k) = sqrt(adjprior.^2 * s2);
   end
else
   E = (cvclass - repmat(id,1,size(cvclass,2))).^2;
   cverr = mean(E,1);
   s2 = sum((E - repmat(cverr,size(E,1),1)).^2,1) / N;
   secverr = sqrt(s2/N);
end

% ----------------------------
function k = findsubtree(Tree,alpha0)
%FINDSUBTREE Find subtree corresponding to specified complexity parameters.

adjfactor = 1 + 100*eps;
alpha = Tree.alpha;
k = zeros(size(alpha0));
for j=1:length(alpha0);
   k(j) = sum(alpha <= alpha0(j)*adjfactor);
end

% -----------------------------
function bestj = selecttree(allalpha,sealpha,treesize)
%SELECTTREE Select the best tree from error rates using some criterion.

% Find the smallest tree that gives roughly the minimum error
[minerr,minloc] = min(allalpha);
if isequal(treesize(1),'m')
   cutoff = minerr * (1 + 100*eps);
else
   cutoff = minerr + sealpha(minloc);
end
j = find(allalpha <= cutoff);
bestj = j(end);

% -----------------------------
function idvec = classname2id(idnames,cnames)
%CLASSNAME2ID Create vector of numeric indices from class name array.

idvec = zeros(length(idnames),1);
if isnumeric(idnames), idnames = cellstr(num2str(idnames)); end
for j=1:length(cnames)
   idvec(strcmp(cnames(j),idnames)) = j;
end

t = find(idvec==0);
if ~isempty(t)
   txt = idnames(t,:);
   if ischar(txt)
      txt = cellstr(txt);
   end
   idvec(t(cellfun('isempty',txt))) = NaN;
end
