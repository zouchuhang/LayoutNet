function [m, maxnum] = mymode(x)
% computes the mode of x (the most commonly occuring value)
% x must be discrete and numeric
% if several values tie for frequency of occurance, each will be returned
% m - the most common value
% maxnum - the number of times that x==m

[g, gn] = grp2idx(x);
maxnum = 0;
m = [];
for i = 1:length(gn)
    count = sum(g==i);
    if count > maxnum
        maxnum = count;
        m = str2num(gn{i});
        c = 1;
    elseif count == maxnum
        m(c+1) = str2num(gn{i});
        c = c + 1;
    end
end

m = sort(m);