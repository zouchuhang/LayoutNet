% crfTrainScript

ncv = 5;
outdir = '/IUS/vmr20/dhoiem/data/mcmcdata';

if 0
    load('/IUS/vmr7/dhoiem/context/data/allimsegs2.mat');
    load('/IUS/vmr7/dhoiem/context/data/rand_indices.mat');
    load([outdir '/traindata5.mat'])
    load('/IUS/vmr7/dhoiem/context/data/mcmcEdgeData.mat');    
end
    


if ~exist('pvSP')
    for tf = 1:numel(cv_images)
        f = cv_images(tf);
        [pvSP{tf}, phSP{tf}, pE{tf}] = mcmcInitialize(spfeatures{f}, efeatures{f}, ...
            adjlist{f}, imsegs(f), vclassifierSP, hclassifierSP, eclassifier, 'labels');        
    end
end


for k = 1:2
    disp(['Iteration: ' num2str(k)]);       
    testind = (floor((k-1)*numel(cv_images)/ncv)+1):(floor(k*numel(cv_images)/ncv));
    trainind = setdiff([1:numel(cv_images)], testind);
    trainind1 = trainind(1:50);                       
    cvtrain1 = cv_images(trainind1);
    if ~exist('vc') || numel(vc)<k
        [tmp1, tmp2, tmp3, tmp4, pgk] = ...
            testMultipleSegmentationsCV2(imsegs(cvtrain1), labdata(trainind1, :), segdata(trainind1, :), smaps(trainind1), ...
            vclassifier(k), hclassifier(k), sclassifier(k), pvSP(trainind1), phSP(trainind1), 1);
        [vc{k}, hc{k}] = splitpg(pgk);
    end
    if ~exist('crfEdgeParam') || numel(crfEdgeParam)<k 
        trainind2 = trainind(51:end); 
        cvtrain2 = cv_images(trainind2);
        crfEdgeParam{k} =  crfTrainInteractionPotential(imsegs(cvtrain2), adjlist(cvtrain2), pE(trainind2));  
        %tmp =  crfTrainInteractionPotential(imsegs(cvtrain2), adjlist(cvtrain2), pE(trainind2)); 
        %crfEdgeParam{k}(:, 3) = tmp(:, 3);
    end
    [crfw{k}, crfpriors{k}] = crfPairwiseTrain3(imsegs(cvtrain1), vc{k}, hc{k}, ...
        adjlist(cvtrain1), pE(trainind1), crfEdgeParam{k}, edgelen(cvtrain1));   
    %crfw4{k} = crfPLTrain2({imsegs(cvtrain1).vert_labels}, {imsegs(cvtrain1).horz_labels}, ...
    %     vc{k}, hc{k}, adjlist(cvtrain1), pE(trainind1), crfEdgeParam{k});    
    disp(num2str([crfw{k}]))
end
save([outdir '/crfTrainParams.mat'], 'vc', 'hc', 'crfEdgeParam', 'crfw', 'crfpriors');


