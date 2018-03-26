function writeAllLabeledImages(imdir, imsegs, pg, outdir)
% writeAllLabeledImages(imdir, imsegs, pg, outdir)
%
% Write label images to specified outdir
%

DO_GT = 0;
if isempty(pg)
    disp('Writing ground truth')
    DO_GT = 1;
end

for f = 1:numel(imsegs)
    
    if DO_GT
        nclasses =  numel(imsegs(f).label_names);
        pg{f} = zeros(imsegs(f).nseg, nclasses);
        lab = imsegs(f).labels;
        ind = find(lab>0);
        pg{f}(ind + imsegs(f).nseg*(lab(ind)-1)) = 1;
    end
    
    im = im2double(imread([imdir '/' imsegs(f).imname]));
    disp([num2str(f) ': ' imsegs(f).imname])
    [vc, hc] = splitpg(pg{f});      
    
    lim = APPgetLabeledImage2(im, imsegs(f), vc, hc);
    
    %system(['cp ' imdir '/' imsegs(f).imname ' ' outdir '/' imsegs(f).imname]);
    imwrite(lim, [outdir '/' strtok(imsegs(f).imname, '.') '.l.jpg'], 'Quality', 90);
end
    
    
    