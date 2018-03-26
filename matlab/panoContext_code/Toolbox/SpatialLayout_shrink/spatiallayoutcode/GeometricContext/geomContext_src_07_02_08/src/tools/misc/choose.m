function z = choose(n, m)
% n choose m = n!/m!/(n-m)!
z = prod([n-m+1:n]) / factorial(m);