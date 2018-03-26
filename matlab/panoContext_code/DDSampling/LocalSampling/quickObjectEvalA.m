function [ scores ] = quickObjectEvalA( objfea, model, typenum )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% objfea = hyps_obj.objfea;
% numobj = length(objfea.sel_hyps);
scores = zeros(length(objfea),typenum);
HYP.sel_hyps = repmat(struct('fixed',false,'selID',[]), typenum, 1);
EMPTY_HYPS = repmat(HYP, length(objfea), 1);
MINSCORE = zeros(length(objfea),10);

model.TRAIN_DATA_IDS = 1:10;
for i = 1:typenum
    ALL_HYPS = EMPTY_HYPS;
    for j = 1:length(objfea)
        ALL_HYPS(j).sel_hyps(i).fixed = true;
        ALL_HYPS(j).sel_hyps(i).selID = j;
    end
    [ sceneImgFea ] = compSceneObjsFeatureA( ALL_HYPS, objfea );
    fea_mask = zeros(1,2870);
    fea_mask((i-1)*224+1:i*224) = 1;
    model.fea_mask = fea_mask;
    [ ~, scores(:,i) ] = sceneHypsEvaluation( sceneImgFea, MINSCORE, model, 11 );
end

end

