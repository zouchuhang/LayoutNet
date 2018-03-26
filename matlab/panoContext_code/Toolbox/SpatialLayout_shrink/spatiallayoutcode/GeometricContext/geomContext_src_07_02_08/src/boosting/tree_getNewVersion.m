function dt_new = tree_getNewVersion(dt)

dt_new=dt;
cut_old = dt.cut;
dt_new.cut=cell(size(cut_old));
dt_new.cut= num2cell(dt.cut);
% dt.cut= num2cell(dt.cut);
 
ind = find(dt.var<0);
for k = 1:numel(ind)
    dt_new.cut{ind(k)} = dt_new.catsplit(cut_old(ind(k)), :);
end
