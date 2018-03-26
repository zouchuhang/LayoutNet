function [ gc_sum, gc_avg, gc_std  ] = gcFeature( gc, indicator, coor2D_indi )
%GCFEATURE Summary of this function goes here
%   Detailed explanation goes here
if size(indicator,1) ~= size(coor2D_indi,1)
    fprintf('Warning: indicator and coor mismatch!\n');
end
% vector_num = size(indicator,1);
data_num = size(indicator, 2);

gc_sum = zeros( data_num, 11);
gc_avg = zeros( data_num, 11);
gc_std = zeros( data_num, 11);

[sH, sW, ~] = size(gc);
% feature = zeros(data_num, 100);

gc = cat(3, gc, sum(gc(:,:,1:4),3), max(gc(:,:,1:4),[],3), sum(gc(:,:,1:6),3), max(gc(:,:,1:6),[],3));
for did = 1:data_num
    vector_valid = indicator(:,did);
    if ~any(vector_valid)
        continue;
    end
    cood_ind = coor2D_indi(vector_valid);    
    loc_orit = zeros(length(cood_ind),11);
    for k = 1:11
        loc_orit(:,k) = gc(cood_ind+(k-1)*sW*sH);
    end
    
    gc_sum(did,:) = sum(loc_orit, 1);
    gc_avg(did,:) = gc_sum(did,:)./length(cood_ind);
    gc_std(did,:) = std(loc_orit, 1);
end

end

