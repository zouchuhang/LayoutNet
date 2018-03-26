% testGraphCutsScript

if 0
    outdir = '../data';
    load([outdir '/rand_indices.mat']);
    load([outdir '/allimsegs2.mat']);
    %load([outdir '/mcmcSuperpixelData.mat']);
    tmp = load('/IUS/vmr20/dhoiem/data/ijcv06/spTrain.mat');
    spfeatures = tmp.spfeatures2;
    vclassifierSP = tmp.spvclassifier2;
    hclassifierSP = tmp.sphclassifier2;
    clear tmp;
    load([outdir '/mcmcSuperpixelClassifier.mat']);
    load([outdir '/mcmcEdgeData.mat']);
    load([outdir '/mcmcEdgeClassifier.mat']);
end

for cvf = 1:numel(cv_images)
    
    disp([num2str(cvf) ': ' imsegs(cv_images(cvf)).imname])
    
    f = cv_images(cvf);    
    c = ceil(cvf/50);
    
    [labv, labh] = testImageGraphCuts(spfeatures{f}, efeatures{f}, adjlist{f}, ...
        vclassifierSP(c), hclassifierSP(c), eclassifier, ecal{c});    

    nsp = numel(labv);
    pg{cvf} = zeros(nsp, 7);
    lab = (labv==1) + (labv==2).*(1+labh) + 7*(labv==3);
    % set most likely label
    pg{cvf}((lab-1)*nsp + [1:nsp]') = 1;
    % set so that h accuracy can be determined
    pg{cvf}((labh)*nsp + [1:nsp]') = 0.5;
end

[vacc, hacc, vcm, hcm] = mcmcProcessResult(imsegs(cv_images), pg);