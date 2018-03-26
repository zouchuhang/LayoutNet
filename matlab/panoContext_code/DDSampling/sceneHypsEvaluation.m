function [ NORM_FEATURE, SCORE ] = sceneHypsEvaluation( FEATURES, MINSCORE, model, aid )
%SCENEHYPSEVALUATION Summary of this function goes here
%   Detailed explanation goes here
TRAIN_DATA_IDS = model.TRAIN_DATA_IDS;
fea_avg = model.fea_avg;
fea_std = model.fea_std;
fea_mask = model.fea_mask;

MINSCORE = min(MINSCORE,4.5);
[B,I] = sort(MINSCORE(:,setdiff(TRAIN_DATA_IDS,aid)),2);
roomfea_sum = cumsum(B(:,1:10),2);
roomfea_prod = cumprod(B(:,1:10),2);
roomfea_algmean = roomfea_sum;
roomfea_geomean = roomfea_prod;
for j = 1:10
    roomfea_algmean(:,j) = roomfea_algmean(:,j)/j;
    roomfea_geomean(:,j) = roomfea_geomean(:,j).^(1/j);
end
ALL_FEATURE = [FEATURES roomfea_sum roomfea_prod roomfea_algmean roomfea_geomean];
fea_num = size(ALL_FEATURE, 1);
NORM_FEATURE = (ALL_FEATURE-repmat(fea_avg, fea_num, 1))./repmat(fea_std, fea_num, 1).*repmat(fea_mask, fea_num, 1);

if ~isempty(model) && isfield(model, 'w')
    w_init = model.w(1:end-1);
    b_init = -model.w(end);

    SCORE = NORM_FEATURE*w_init' - b_init;
else
    SCORE = [];
end

end

