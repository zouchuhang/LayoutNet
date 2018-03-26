function [blocks,dag] = graphComponents(A)
% GRAPHCOMPONENTS : Connected or strongly connected components of a graph.
%
% blocks = graphComponents(A);
% [blocks,dag] = graphComponents(A);
%
% Input:  A         is the n by n adjacency matrix for an n-vertex graph.
%
% Output: blocks    is an n-vector of integers 1:k, where A has k components, 
%                   that labels the vertices of A according to component.
%         dag       (optional) is the k-vertex acyclic directed graph obtained 
%                   by contracting each component of A into a vertex,
%                   with component sizes on the diagonal and contracted-edge
%                   counts off the diagonal.
%
% If the input A is undirected (i.e. symmetric), the blocks are its connected
% components and the dag has no edges (i.e. is a diagonal matrix).
%
% If the input A is directed (i.e. unsymmetric), the blocks are its strongly 
% connected components, numbered in topological order according to the dag.
%
% See also DMPERM.
%
% John Gilbert, Xerox PARC, 8 June 1992.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% http://www.cerfacs.fr/algor/Softs/MESHPART/

 
[n,m] = size(A);
if n ~= m, error ('Adjacency matrix must be square'), end;

if ~all(diag(A)) 
    [p,p,r,r] = dmperm(A|speye(size(A)));
else
    [p,p,r,r] = dmperm(A);  
end;
 
% Now the i-th component of A(p,p) is r(i):r(i+1)-1.

sizes = diff(r);        % Sizes of components, in vertices.
k = length(sizes);      % Number of components.
 
% Now compute an array "blocks" that maps vertices of A to components;
% First, it will map vertices of A(p,p) to components...
 
blocks = zeros(1,n);
blocks(r(1:k)) = ones(1,k);
blocks = cumsum(blocks);
 
% Second, permute it so it maps vertices of A to components.
 
blocks(p) = blocks;
 
% Compute the acyclic condensation by using the map,
% taking care to get the diagonal right.
 
if nargout > 1
    [i,j] = find(A);
    dag = sparse(blocks(i),blocks(j),1,k,k);
    dag = dag + diag(sparse(sizes)) - diag(diag(dag));
end;
