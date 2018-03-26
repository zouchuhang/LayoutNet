function ordered_textons = select_diverse_textons(tsim, nt)

next_t = zeros(nt, 1);
vals = zeros(nt, 1);
indices = (1:size(tsim, 1));
for t = 1:nt
    total_dsim = mean(tsim, 1);
    [vals(t), index] = max(total_dsim);
    next_t(t) = indices(index);
    tsim = tsim([(1:index-1) (index+1:end)], [(1:index-1) (index+1:end)]);
    indices = indices([(1:index-1) (index+1:end)]);
end

ordered_textons = next_t;