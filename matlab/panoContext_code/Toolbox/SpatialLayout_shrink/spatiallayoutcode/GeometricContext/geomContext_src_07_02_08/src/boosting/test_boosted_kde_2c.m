function p = test_boosted_kde_2c(density, x, data)
% used to evaluate the likelihood of a point from a kernel density estimate
% density is the density function
% x are the points at which f is defined (assumed to be equally spaced)
% data are the data points to be evaluated
% p is the value of f(xi) where x(xi) is the closest value in x to y

n = length(x);
wx = x(2)-x(1);
indices = round((data- x(1))/wx)+1;
indices = min(max(indices, 1), n);
p = density(indices);