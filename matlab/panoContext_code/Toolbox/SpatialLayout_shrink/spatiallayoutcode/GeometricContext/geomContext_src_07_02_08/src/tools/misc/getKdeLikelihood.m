function p = getKdeLikelihood(f, x, y)
% p = getKdeLikelihood(f, x, y)
% Evaluates the likelihood of a point from a kernel density estimate
%
% Input:
%   f: the density function
%   x: the points at which f is defined (assumed to be equally spaced)
%   y: the data points to be evaluated
% Output:
%   p: the value of f(xi) where x(xi) is the closest value in x to y

n = length(x);
wx = x(2)-x(1);
indices = round((y - x(1))/wx)+1;
indices = min(max(indices, 1), n);
p = f(indices);