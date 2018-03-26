function [nnidx, dists] = annsearch(X0, X, k, varargin)
%ANNSEARCH Approximate Nearest Neighbor Search
%
% $ Syntax $
%   - nnidx = annsearch(X0, [], k)
%   - nnidx = annsearch(X0, X, k)
%   - nnidx = annsearch(X0, [], k, ...)
%   - nnidx = annsearch(X0, X, k, ...)
%   - [nnidx, dists] = annsearch(...)
%
% $ Arguments $
%   - X0:       the matrix of points for constructing the KD tree
%   - X:        the matrix of points for query
%   - k:        the number of neighbors for each query point
%   - nnidx:    the indices of neighboring points
%   - dists:    the distances between the query points to their neighbors
%
% $ Description $
%   - nnidx = annsearch(X0, [], k) approximately searches the neighboring
%     points of each point of X0 within X0. In the findings, the query
%     point itself is excluded from its neighbor set. k non-trivial 
%     neighbors for each query point is searched.
%   
%   - nnidx = annsearch(X0, X, k) approximately searches the neighbors of
%     points in X within the point set specified by X0. Note that even
%     though there are points in X0 which are exactly some point in X, they
%     will not be excluded from the results.
%
%   - nnidx = annsearch(X0, [], k, ...)
%   - nnidx = annsearch(X0, X, k, ...) conducts the search with
%     user-customized options. The options that can be set are given
%     as follows:
%     \*
%     \t  Table 1. Options for ANN search
%     \h     name       &      description      \\
%           'errbound'  &  The allowable upper bound on error (default = 0) \\
%           'split'     &  The rule of splitting (default = 'suggest') \\
%           'search'    &  The search method (default = 'normal') \\
%     \*
%    There are following splitting rule in constructing the KD tree:
%     \*
%     \t  Table 2. The Splitting rules for KD tree construction
%     \h     name        &   description   \\
%            'std'       &  the standard optimized kd-splitting rule \\
%            'midpt'     &  the midpoint splitting rule \\
%            'fair'      &  the fair splitting rule \\
%            'sl_midpt'  &  the sliding midpoint splitting rule \\
%            'sl_fair'   &  the sliding fair splitting rule \\
%            'suggest'   &  the author's suggestion for best (default) \\
%     \*
%   There are following search methods:
%     \*
%     \t  Table 3. The search methods
%     \h     name       &     description  \\
%            'normal'   &    the normal search method (default) \\
%            'priority' &    the priority search method
%     \*
%
%   - [nnidx, dists] = annsearch(...) also returns the distances between
%     the query points and their neighbors.
%
% $ Remarks $
%   - In X0 or X, each column represents a sample point.
%   - If there are n points in the query set, then both nnidx and dists
%     would be k x n matrix, with each column recording the neighbor
%     information about the corresponding point.
%
% $ History $
%   - Created by Dahua Lin on Apr 21, 2006
%

%% parse and verify arguments

% for argument number
if nargin < 3
    error('annerror:invalidarg', ...
        'The number of input arguments should not be less than 3 for annsearch');
end
if nargout == 0
    return;
elseif nargout > 2
    error('annerror:invalidarg', ...
        'The number of output arguments should not be larger than 2 for annsearch');
end

% for base-target relation
if ~isempty(X)
    exclude_self = false;
else
    X = X0;
    k = k + 1;
    exclude_self = true;
end

% for sizes
[d, n0] = size(X0);
if size(X, 1) ~= d
    error('annerror:invalidarg', ...
        'The dimension of training and query points should be consistent');
end
if k >= n0
    error('annerror:invalidarg', ...
        'The value k (neighborhood size) should be less than n0, the size of the whole set');
end

% for options
opts.errbound = 0;
opts.split = 'suggest';
opts.search = 'normal';
opts = slparseprops(opts, varargin{:});



%% invoke the wrapped core function

if nargout == 1

    nnidx = annsearch_wrapper( ...
        X0, ...
        X, ...
        k, ...
        opts.errbound, ...
        get_splitrule_id(opts.split), ...
        get_searchmethod_id(opts.search));

else
    
    [nnidx, dists] = annsearch_wrapper( ...
        X0, ...
        X, ...
        k, ...
        opts.errbound, ...
        get_splitrule_id(opts.split), ...
        get_searchmethod_id(opts.search));

end

%% post-processing

if exclude_self
    nnidx = nnidx(2:end, :);
    
    if nargout >= 2
        dists = dists(2:end, :);
    end
end

nnidx = nnidx + 1;  % from zero-based to one-based

if nargout >= 2     % from square distance to euclidean distance
    dists = sqrt(dists);
end


%% Auxiliary functions

function id = get_splitrule_id(s)

switch s
    case 'suggest'
        id = 5;
    case 'std'
        id = 0;
    case 'midpt'
        id = 1;
    case 'fair'
        id = 2;
    case 'sl_midpt'
        id = 3;
    case 'sl_fair'
        id = 4;
    otherwise
        error('annerror:invalidarg', ...
            'Invalid split rule %s for annsearch', s);
end


function id = get_searchmethod_id(s)

switch s
    case 'normal'
        id = 0;
    case 'priority'
        id = 1;
    otherwise
        error('annerror:invalidarg', ...
            'Invalid search method %s for annsearch', s);
end







