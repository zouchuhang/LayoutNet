function y = randsample(n, k)
%RANDSAMPLE Random sampling, without replacement
%   Y = RANDSAMPLE(N,K) returns K values sampled at random, without
%   replacement, from the integers 1:N.

%   Copyright 1993-2002 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2002/03/13 23:15:54 $

% RANDSAMPLE does not (yet) implement weighted sampling.

if nargin < 2
    error('Requires two input arguments.');
end

% If the sample is a sizeable fraction of the population, just
% randomize the whole population (which involves a full sort
% of n random values), and take the first k.
if 4*k > n
    rp = randperm(n);
    y = rp(1:k);

% If the sample is a small fraction of the population, a full
% sort is wasteful.  Repeatedly sample with replacement until
% there are k unique values.
else
    x = zeros(1,n); % flags
    sumx = 0;
    while sumx < k
        x(ceil(n * rand(1,k-sumx))) = 1; % sample w/replacement
        sumx = sum(x); % count how many unique elements so far
    end
    y = find(x > 0);
    y = y(randperm(k));
end

% a scalar loop version
%
% x = 1:n;
% n = n:(-1):(n-k+1);
% y = zeros(1,k);
% j = ceil(n .* rand(1,k));
% for i = 1:k
%     y(i) = x(j(i));
%     x(j(i)) = x(n(i));
% end
