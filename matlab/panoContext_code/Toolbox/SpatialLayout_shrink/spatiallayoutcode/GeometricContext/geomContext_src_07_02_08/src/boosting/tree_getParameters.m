function [var, cut, children, catsplit] = tree_getParameters(dt)

if iscell(dt.cut)
    var = dt.var;    
    children = dt.children;
    
    cut = zeros(size(dt.cut));
    catind = var<0;
    ncat = sum(catind);  
    catsplit = cell(ncat, 2);
    ncat = 0;
    for k = 1:numel(cut)
        if var(k)>=0
            cut(k) = dt.cut{k};
        else
            ncat = ncat+1;
            cut(k) = ncat;
            catsplit(ncat, :) = dt.cut{k};
        end
    end      
else
    var = dt.var;
    cut = dt.cut;
    children = dt.children;
%     catsplit = dt.catsplit;
if any(var<0)
        catsplit = dt.catsplit;
    else
        catsplit = cell(0,2);
    end
end
