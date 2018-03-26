function [ cn_hist, cn_entropy ] = colorNameFeature( cn, indicator, coor2D_indi )
%COLORNAMEFEATURE Summary of this function goes here
%   Detailed explanation goes here
if size(indicator,1) ~= size(coor2D_indi,1)
    fprintf('Warning: indicator and coor mismatch!\n');
end

data_num = size(indicator, 2);
cn_hist = zeros(data_num, 11);
% cn_entropy = zeros(data_num,1);
for did = 1:data_num
    vector_valid = indicator(:,did);
    if ~any(vector_valid)
        continue;
    end
    cood_ind = coor2D_indi(vector_valid);    
    loc_cn = cn(cood_ind);
    cn_hist(did,:) = hist(loc_cn, 1:11)./length(loc_cn);       
end

cn_entropy = -sum(cn_hist.*log2(cn_hist+0.0001), 2);

end

