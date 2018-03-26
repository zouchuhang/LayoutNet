function pval=statctexact(x,wts,tstar,dispopt)
%STATCTEXACT Compute exact p-value for contingency table
%   P=STATCTEXACT(X,WTS,T,DISPOPT) uses a network algorithm to compute
%   the exact p-value P for a 2-by-K contingency table X.  The test 
%   statistic T is the weighted sum of the elements in the first row.
%   Set DISPOPT=true to display debugging output.
%
%   Private function used by the RANKSUM function.

%   Copyright 2003-2004 The MathWorks, Inc. 
%   $Revision: 1.1.6.2 $  $Date: 2004/01/24 09:36:26 $


[r,c] = size(x);
if (r~=2)
   error('Internal error, table must have two rows.');
end

if (nargin<2), wts = []; end
if (nargin<3), tstar = []; end
if (nargin<4), dispopt = false; end

% Get the test statistic and weights if not already done
if isempty(wts) || isempty(tstar)
   [tstar, wts, expected] = teststat(x, wts, dispopt);
else
   expected = [];   % will not be used
end

% Start of network algorithm for exact p-value computation

% Make nodes and arcs
[nodes,arcs] = makenodes(x,wts,expected);

% Backward induction
nodes = backward(nodes,arcs);

% Forward pass
pvals = getsigprob(nodes, arcs, tstar, dispopt);

TP = nodes{4,1};
pvals = pvals / TP;
pval = pvals(2) + min(pvals(1), pvals(3));
p2 = min(1, 2 * pval);

if dispopt
   disp(        'Exact results:');
   disp(sprintf('   test statistic = %g', tstar));
   disp(sprintf('   Prob[T <  %g] = %g', tstar, pvals(1)));
   disp(sprintf('   Prob[T <= %g] = %g', tstar, pvals(1)+pvals(2)));
   disp(sprintf('   Prob[T =  %g] = %g', tstar, pvals(2)));
   disp(sprintf('   Prob[T >= %g] = %g', tstar, pvals(2)+pvals(3)));
   disp(sprintf('   Prob[T >  %g] = %g', tstar, pvals(3)));
   disp(sprintf('   2-sided p = %g', p2));
end

% ----------------------------------------------------------
function [tstar, wts, expected] = teststat(x, wts, dispopt);
% Compute the test statistic for this observed table

[r,c] = size(x);
rowsum = sum(x,2);
colsum = sum(x,1);
if (length(wts) == 0)
   obs = x;
   expected = repmat(rowsum,1,c) .* repmat(colsum,r,1) ./ sum(rowsum);
   tstar = sum(sum((obs-expected).^2 ./ expected));
else
   expected = [];
   tstar = sum(wts .* x(1,:));
end
wts = wts(:)';

% ---------------------------------------------------
function [nodes,arcs] = makenodes(x,wts,expected)
%MAKENODES  Make structures describing nodes and arcs
%    [nodes,arcs] = makenodes(x,wts,expected)

[r,c] = size(x);
rowsum = sum(x,2);
colsum = sum(x,1);
oldnodes = zeros(1,2); % nodes added during the last pass
oldlo = 0;             % min possible sum so far
oldhi = 0;             % max possible sum so far
oldnn = 1;             % node numbers (row numbers) from last pass
xsum = rowsum(1);      % sum of entries in first row
nodecount = 1;

nodes = cell(4,c+1);   % to hold nodes
%      row 1:  n-by-2 array, n = # of nodes, row = [j,mj]
%      row 2:  n-vector of longest path to end from here
%      row 3:  n-vector of shortest path to end from here
%      row 4:  n-vector of total probability to end from here

arcs = cell(3,c);      % to hold node connections (arcs) in the network
%      row 1:  n-by-2 array, n = # of connections, row = pair connected
%      row 2:  n-vector of arc lengths
%      row 3:  n-vector of arc probabilities

nodes{1,1} = zeros(1,2);
nodes{2,c+1} = 0;
nodes{3,c+1} = 0;
nodes{4,c+1} = 1;
for j=1:c            % loop over nodes
   % Figure out which nodes are possible at the next step
   nj = colsum(j);
   lo = max(oldlo, xsum-sum(colsum(j+1:end)));
   hi = min(xsum, oldhi+nj);
   newnodes = zeros(hi-lo+1,2);
   newnodes(:,1) = j;
   newnodes(:,2) = (lo:hi)';
   newnn = 1:size(newnodes,1);
   nodecount = nodecount + size(newnodes,1);
   nodes{1,j+1} = newnodes;
   
   % Figure out which arcs are possible to the next step
   [a0,a1] = meshgrid(oldnn,newnn);
   a0 = a0(:);
   a1 = a1(:);
   oldsum = oldnodes(a0,2);
   newsum = newnodes(a1,2);
   xj = newsum - oldsum;
   ok = (xj >= 0) & (xj <= nj);
   arcs{1,j} = [a0(ok) a1(ok)];  % arc connections
   xj = xj(ok);
   if (length(wts) > 0)          % arch lengths
      arcs{2,j} = wts(j)*xj;
   else
      arcs{2,j} = 2 * (xj - expected(1,j)).^2 ./ expected(1,j);
   end
   pj = exp(gammaln(nj+1) - gammaln(xj+1) - gammaln(nj-xj+1));
   arcs{3,j} = pj;               % arc probabilities
   
   % Update data structures
   oldlo = lo;
   oldhi = hi;
   oldnodes = newnodes;
   oldnn = newnn;
end

% -----------------------------------------------------------
function nodes = backward(nodes,arcs)
%BACKWARD  Do backward induction, add information to NODES array
%   nodes = backward(nodes,arcs)

% initialize for final node
c = size(nodes,2) - 1;
startSP = zeros(1);
startLP = startSP;
startTP = ones(1);
startnode = nodes{1,c+1};
for j=c:-1:1
   % destination nodes are previous start nodes
   endSP = startSP;
   endLP = startLP;
   endTP = startTP;
   endnode = startnode;
   
   % get new start nodes and information about them
   a = arcs{1,j};
   startnode = nodes{1,j};
   startmax = max(a(:,1));
   startSP = zeros(startmax,1);
   startLP = startSP;
   startTP = startSP;
   arclen = arcs{2,j};
   arcprob = arcs{3,j};
   for nodenum=1:startmax
      % for each start node, compute SP, LP, TP
      k1 = find(a(:,1) == nodenum);
      k2 = a(k1,2);
      startLP(nodenum) = max(arclen(k1) + endLP(k2));
      startSP(nodenum) = min(arclen(k1) + endSP(k2));
      startTP(nodenum) = sum(arcprob(k1) .* endTP(k2));
   end
   
   % store information about nodes at this level
   nodes{2,j} = startLP;
   nodes{3,j} = startSP;
   nodes{4,j} = startTP;
end

% ----------------------------------------------------
function pvals = getsigprob(nodes, arcs, tstar, dispopt)
%GETSIGPROB Get p-values by scanning the network

NROWS = 50;

pvals = zeros(3,1);    % [Prob<T, Prob=T, Prob>T]
stack = zeros(NROWS, 4);
stack(:,1) = Inf;
stack(1,1) = 1;        % level of current node
stack(1,2) = 1;        % number at this level of current node
stack(1,3) = 0;        % length so far to this node
stack(1,4) = 1;        % probability so far of reaching this node
N = size(stack,1);

i1 = 0; i2 = 0; i3 = 0;

while(1)
   % Get next node to process, visiting lowest levels first
   minlevel = min(stack((stack(1:N)>0)));
   if (isinf(minlevel)), break; end
   sp = find(stack(1:N)==minlevel);
   sp = sp(1);
   
   L = stack(sp,1);
   J = stack(sp,2);
   pastL = stack(sp,3);
   pastP = stack(sp,4);
   stack(sp,1) = Inf;
   
   % Get info for arcs at level L and their target nodes
   nj = nodes{1,L+1};
   LP = nodes{2,L+1};
   SP = nodes{3,L+1};
   TP = nodes{4,L+1};
   aj = arcs{1,L};
   arclen = arcs{2,L};
   arcprob = arcs{3,L};

   % Look only at arcs from node J
   seps = sqrt(eps);
   arows = find(aj(:,1)==J)';
   for k=arows
      tonode = aj(k,2);
      thisL = arclen(k);
      thisP = pastP * arcprob(k);
      len = pastL + thisL;

      % See if no paths from here are signicant
      if (len + LP(tonode) < tstar - seps)
         pvals(1) = pvals(1) + thisP * TP(tonode);
         continue;

      % See if all paths from here are significant
      elseif (len + SP(tonode) > tstar + seps)
         pvals(3) = pvals(3) + thisP * TP(tonode);
         continue;

      % See if there is no range, then we match exactly
      elseif (SP(tonode) == LP(tonode))
         pvals(2) = pvals(2) + thisP * TP(tonode);
         continue;

      % See if we can merge this with another already stored
      else
         % Find a stored node that matches this one
         r = find(stack(:,1) == L+1);
         if (any(r))
            r = r(stack(r,2) == tonode);
            if (any(r))
               r = r(abs(stack(r,3) - len) < seps);
            end
         end
         
         if (any(r))
            % If one is found, merge this one with it
            sp = r(1);
            stack(sp,4) = stack(sp,4) + thisP;
            i1 = i1+1;
         else
            % Otherwise add a new node, extending array if necessary
            z = find(isinf(stack(:,1)));
            if (isempty(z))
               i2 = i2+1;
               block = zeros(NROWS,4);
               block(:,1) = Inf;
               stack = [stack; block];
               sp = N+1;
               N = N+NROWS;
            else
               i3 = i3+1;
               sp = z(1);
            end
            stack(sp,1) = L+1;
            stack(sp,2) = tonode;
            stack(sp,3) = len;
            stack(sp,4) = thisP;
         end
      end
   end
end

if dispopt
   disp('merged, extended, inserted = ');
   disp([i1 i2 i3]);
end
