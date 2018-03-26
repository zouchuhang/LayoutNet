% testGraphCutsScript3

if 0
    imdir = '../images/all_images';
    outdir = '../data';
    load([outdir '/rand_indices.mat']);
    load([outdir '/allimsegs2.mat']);
    load([outdir '/mcmcSuperpixelClassifier.mat']);
    load([outdir '/mcmcEdgeData.mat']);
    load([outdir '/mcmcEdgeClassifier.mat']);
end

alpha = [1:1:3];

for k = 1:numel(alpha);

    disp(num2str(alpha(k)))
    
    for cvf = 1:numel(cv_images)

        disp([num2str(cvf) ': ' imsegs(cv_images(cvf)).imname])        
        
        f = cv_images(cvf);    
        c = ceil(cvf/50);

        im = im2double(imread([imdir '/' imsegs(f).imname]));
        wseg = segmentWatershed(im, 1);
        
        cimages = pg2confidenceImages(imsegs(f), pg(cvf));
        [labv, labh] = testImageGraphCuts3(im, wseg, cimages{1}, alpha(k));    

        nsp = numel(labv);
        pg2{cvf} = zeros(nsp, 7);
        lab = (labv==1) + (labv==2).*(1+labh) + 7*(labv==3);
        % set most likely label
        pg2{cvf}((lab-1)*nsp + [1:nsp]') = 1;
        % set so that h accuracy can be determined
        pg2{cvf}((labh)*nsp + [1:nsp]') = 0.5;
        
        [pv1, ph1] = splitpg(pg{cvf});
        lim1 = APPgetLabeledImage2(im, imsegs(f), pv1, ph1);
        [pv2, ph2] = splitpg(pg2{cvf});
        lim2 = APPgetLabeledImage2(im, wseg, pv2, ph2);        
        figure(1), hold off, imagesc(lim1), axis image;
        figure(2), hold off, imagesc(lim2), axis image;
        drawnow; pause(1)
    end

    [vaccgc(k), haccgc(k)] = mcmcProcessResult(imsegs(cv_images), pg2);
    disp(num2str([vaccgc ; haccgc]))
end