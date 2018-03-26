function [pg, energy, map, allenergy] = mcmcCvEvaluateSA(imsegs, imdir, vclassifierSP, hclassifierSP, ...
    eclassifier, vclassifier, hclassifier, segclassifier, spdata, adjlist, edata)

ncv = numel(vclassifier);   

niter = 5;

vacc = zeros(numel(imsegs), 1);
hacc = zeros(numel(imsegs), 1);

%load '../data/tmp.mat'

firstf = 1;
% disp('Temporarily starting from 51!')
% firstf = 51;
% load '../data/tmp.mat'
% size(pg)

for f = firstf:numel(imsegs)
   
  
    
    k = ceil(f/numel(imsegs)*ncv);
    
    disp([num2str(f) ': ' imsegs(f).imname])
    
    im = im2double(imread([imdir '/' imsegs(f).imname]));
    
    [pg{f}, map{f}, allenergy{f}] = mcmcTestImageSA(im, imsegs(f), ...
        vclassifierSP, hclassifierSP, eclassifier, vclassifier(k), hclassifier(k), ...
        segclassifier(k), niter, spdata{f}, adjlist{f}, edata{f});
                 
   [vacc(f), hacc(f), vtotal(f), htotal(f)] = mcmcProcessResult(imsegs(f), pg(f));
   energy(f) = allenergy{f}(end);     
   
    %save '../data/tmp.mat' pg map vacc hacc vtotal htotal 

    %disp('***** Report *****') ;
    disp('Total')
    disp(['vacc: ' num2str(sum(vacc)/sum(vtotal)) '   hacc: ' num2str(sum(hacc)/sum(htotal))])
    %disp('*** End Report ***');
    
end





