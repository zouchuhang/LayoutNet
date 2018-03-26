function imsegs2 = APPtransferLabels(imsegs1, imsegs2, l2v, l2h)
%
% imsegs2 = transferLabels(imsegs1, imsegs2)
%
% Transfers labeling from imsegs1 to imsegs2.  imsegs1.segimage and ground
% truth labels are used to create a label map, which is then used to label
% the superpixels in imsegs2.  l2v and l2h give conversions from full label
% index to vertical and horz sublabel indices, respectively.
% (default: l2v = [1 2 2 2 2 2 3], l2h = [0 1 2 3 4 5 0])

if ~exist('l2v', 'var')
    l2v = [1 2 2 2 2 2 3];
end
if ~exist('l2h', 'var')
    l2h = [0 1 2 3 4 5 0];
end

for f = 1:numel(imsegs1)
    
    disp(num2str(f))
    
    stats = regionprops(imsegs1(f).segimage, 'PixelIdxList');
    idx = {stats.PixelIdxList};
    lmap = zeros(size(imsegs1(f).segimage));
    for s = 1:imsegs1(f).nseg
        lmap(idx{s}) = imsegs1(f).labels(s);
    end
    
    % set imsegs2
    nseg = imsegs2(f).nseg;
    imsegs2(f).vlabels = cell(nseg, 1);
    imsegs2(f).hlabels = cell(nseg, 1);
    imsegs2(f).labels = zeros(nseg, 1);
    imsegs2(f).vert_labels = zeros(nseg, 1);
    imsegs2(f).horz_labels = zeros(nseg, 1);
    imsegs2(f).label_names = imsegs1(f).label_names;
    imsegs2(f).vert_names = imsegs1(f).vert_names;
    imsegs2(f).horz_names = imsegs1(f).horz_names;
    
    stats = regionprops(imsegs2(f).segimage, 'PixelIdxList');
    idx = {stats.PixelIdxList};
        
    for s = 1:nseg
        mcl = mode(lmap(idx{s}));
        
        imsegs2(f).labels(s) = mcl;
        if mcl==0     
            imsegs2(f).vert_labels(s) = 0;
            imsegs2(f).horz_labels(s) = 0;
        else            
            imsegs2(f).vert_labels(s) = l2v(mcl);
            imsegs2(f).horz_labels(s) = l2h(mcl);
        end
        
        if imsegs2(f).vert_labels(s)==0
            imsegs2(f).vlabels{s} = '---';
        else
            imsegs2(f).vlabels{s} = imsegs2(f).vert_names{imsegs2(f).vert_labels(s)};
        end
        if imsegs2(f).horz_labels(s)==0
            imsegs2(f).hlabels{s} = '---';
        else
            imsegs2(f).hlabels{s} = imsegs2(f).horz_names{imsegs2(f).horz_labels(s)};
        end           
    end
end
figure(1), imagesc(imsegs1.vert_labels(imsegs1.segimage)), axis image, colormap jet        
        