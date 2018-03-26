% testGraphCutsScript2

if 0
    outdir = '../data';
    load([outdir '/rand_indices.mat']);
    load([outdir '/allimsegs2.mat']);
    load([outdir '/mcmcSuperpixelClassifier.mat']);
    load([outdir '/mcmcEdgeData.mat']);
    load([outdir '/mcmcEdgeClassifier.mat']);
end

alpha = [0:0.5:3];

for k = 1:numel(alpha);

    disp(num2str(alpha(k)))
    
    for cvf = 1:numel(cv_images)

        %disp([num2str(cvf) ': ' imsegs(cv_images(cvf)).imname])

        f = cv_images(cvf);    
        c = ceil(cvf/50);

        [labv, labh] = testImageGraphCuts2(pg{cvf}, efeatures{f}, adjlist{f}, eclassifier, ecal{c}, alpha(k));    

        nsp = numel(labv);
        pg2{cvf} = zeros(nsp, 7);
        lab = (labv==1) + (labv==2).*(1+labh) + 7*(labv==3);
        % set most likely label
        pg2{cvf}((lab-1)*nsp + [1:nsp]') = 1;
        % set so that h accuracy can be determined
        pg2{cvf}((labh)*nsp + [1:nsp]') = 0.5;
    end

    [vaccgc(k), haccgc(k)] = mcmcProcessResult(imsegs(cv_images), pg2);
    disp(num2str([vaccgc ; haccgc]))
end