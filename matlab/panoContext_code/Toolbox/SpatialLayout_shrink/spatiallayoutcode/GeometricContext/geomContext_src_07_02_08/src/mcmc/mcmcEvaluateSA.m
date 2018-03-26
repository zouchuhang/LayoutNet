function [pg, energy, map, allenergy] = mcmcEvaluate(imsegs, imdir, vclassifierSP, hclassifierSP, ...
    eclassifier, ecal, vclassifier, hclassifier, sclassifier, spdata, adjlist, edata)

ncv = numel(vclassifier);

niter = 10;

vacc = zeros(numel(imsegs), 1);
hacc = zeros(numel(imsegs), 1);

%load '../data/tmp.mat'
  
for f = 1:numel(imsegs)
    
    c = ceil(f/numel(imsegs)*ncv);
    
    disp([num2str(f) ': ' imsegs(f).imname])
    
    im = im2double(imread([imdir '/' imsegs(f).imname]));
    
    [pg{f}, map{f}, allenergy{f}] = mcmcTestImageSA(im, imsegs(f), ...
        vclassifierSP, hclassifierSP, eclassifier, ecal, vclassifier(c), hclassifier(c), ...
        sclassifier(c), niter, spdata{f}, adjlist{f}, edata{f});
                    
    %energy(f) = allenergy{f}(end);     
   
    %save '../data/tmp.mat' pg vacc hacc vtotal htotal energy

    %disp('***** Report *****') ;
    %disp('Total')
    %disp(['vacc: ' num2str(sum(vacc)/sum(vtotal)) '   hacc: ' num2str(sum(hacc)/sum(htotal))])
    %disp('*** End Report ***');
    
end

[vacc(f), hacc(f), vtotal(f), htotal(f)] = mcmcProcessResult(imsegs(f), pg(f));



