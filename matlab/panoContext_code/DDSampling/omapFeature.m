function [ om_sum, om_avg, om_std ] = omapFeature( omap, indicator, coor2D_indi )
%OMAPFEATURE Summary of this function goes here
%   Detailed explanation goes here
if size(indicator,1) ~= size(coor2D_indi,1)
    fprintf('Warning: indicator and coor mismatch!\n');
end
% vector_num = size(indicator,1);
data_num = size(indicator, 2);

om_sum = zeros( data_num, 5);
om_avg = zeros( data_num, 5);
om_std = zeros( data_num, 5);

[sH, sW, ~] = size(omap);
% feature = zeros(data_num, 100);

omap = cat(3, omap, sum(omap(:,:,2:3),3), max(omap(:,:,2), omap(:,:,3)));
for did = 1:data_num
    vector_valid = indicator(:,did);
    if ~any(vector_valid)
        continue;
    end
    cood_ind = coor2D_indi(vector_valid);    
    loc_orit = zeros(length(cood_ind),5);
    for k = 1:5
        loc_orit(:,k) = omap(cood_ind+(k-1)*sW*sH);
    end
    
    om_sum(did,:) = sum(loc_orit, 1);
    om_avg(did,:) = om_sum(did,:)./length(cood_ind);
    om_std(did,:) = std(loc_orit, 1);
end

end

