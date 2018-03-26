function [pg, pgIter] = mcmcEvaluate(imsegs, imdir, vclassifierSP, hclassifierSP, ...
    eclassifier, vclassifier, hclassifier, segclassifier, spdata, adjlist, edata)

niter = 3000;

istep = 50;

vacc = zeros(numel(imsegs), floor(niter/istep)+1);
hacc = zeros(numel(imsegs), floor(niter/istep)+1);

%load '../data/tmp.mat'
  
for f = 1:numel(imsegs)
    
    disp([num2str(f) ': ' imsegs(f).imname])
    
    im = im2double(imread([imdir '/' imsegs(f).imname]));
    
    [pg{f}, pgIter{f}] = mcmcTestImage(im, imsegs(f), ...
        vclassifierSP, hclassifierSP, eclassifier, vclassifier, hclassifier, ...
        segclassifier, niter, spdata{f}, adjlist{f}, edata{f});
    %pgIter{f} = pgIter{f}(:, :, [1 istep:istep:end]);                    
    c = 0;
    for k = [1 istep:istep:niter]
                
        c = c + 1;
        [vacc(f,c), hacc(f,c), vtotal(f), htotal(f)] = ...
            mcmcProcessResult(imsegs(f), pgIter(f), k);
        
        
%         tmppg = mean(pgIter{f}(:, :, 1:k), 3);
%     
%         [maxval, maxlab] = max(tmppg, [], 2);
%         maxlab = maxlab(imsegs(f).segimage);
%         labels = imsegs(f).labels(imsegs(f).segimage);
%         vacc = sum(sum((maxlab==1 & labels==1) | (maxlab==7 & labels==7) ...
%             | (maxlab>1 & maxlab < 7 & labels>1 & labels<7)));
%     
%         [maxval, maxlab] = max(tmppg(:, 2:6), [], 2);
%         maxlab = maxlab(imsegs(f).segimage);
%         labels = imsegs(f).horz_labels(imsegs(f).segimage);
%         hacc(c,f) = sum(sum((maxlab==labels))) / sum(labels(:)~=0);
        
    end
    pgIter{f} = pgIter{f}(:, :, [1 istep:istep:end]);    
    save '../data/tmp.mat' pgIter pg vacc hacc
    figure(4), hold off, plot([1 istep:istep:niter], sum(vacc(1:f, :), 1)/sum(vtotal(1:f))), axis([1 niter 0.4 0.95])
    figure(5), hold off, plot([1 istep:istep:niter], sum(hacc(1:f, :), 1)/sum(htotal(1:f))), axis([1 niter 0.3 0.8])
    pause(1);
    
end





