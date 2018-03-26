function imsegs = labelImageGeometry(im, imsegs)
% Label the superpixels of the image

segimage = double(imsegs.segimage);
nc = 7;

if max(size(segimage)>600)
    rs = 600/max(size(segimage));
    segimage = imresize(segimage, rs, 'nearest');
    im = imresize(im, rs, 'nearest');
end

try
    labels = imsegs.labels;
catch
    labels = zeros(imsegs.nseg, 1);
end

[imh, imw] = size(segimage);
grayim = repmat(rgb2gray(im), [1 1 3]);

stats = regionprops(segimage, 'Centroid');

yim = repmat((1:imh)', 1, imw);
xim = repmat((1:imw), imh, 1);

centroids = reshape([stats(:).Centroid], [2 imsegs.nseg])';
meanx = centroids(:, 1);
meany = centroids(:, 2);

nseg = imsegs.nseg;

%mapim = label2rgb([1:7]', 'jet', [1 1 1]);
mapim = hsv2rgb([1:7]'/7, ones(7, 1), ones(7, 1));
figure(2), imagesc(mapim)

figure(3), imagesc(im), axis image

tmplim = labels(segimage); 
hue = tmplim / nc;
lim = hsv2rgb(hue, hue>0, ones(size(hue)))*0.5 + 0.5*grayim;
%tmplim(1:7) = [1:7];
%lim = im2double(label2rgb(tmplim, 'jet', [1 1 1]))*0.5 + repmat(grayim, [1 1 3])*0.5;
figure(1), imagesc(lim), axis image


x = 1;
while (1)  
        
    %disp(['num unlabeled = ' num2str(sum(labels==0))])
    [x, y, b] = ginput(1);
    x = round(x);
    y = round(y);
        
    if b=='q'
        break;
    end         

    if b=='z'            
        [tmp1, tmp2, b] = ginput(1);
        lab = str2num(char(b));
        if lab>0 && lab<=nc
            ind = find(labels==0);
            labels(ind) = lab;
            break;
        end                
    end

    tind = find(y>=1 & y<=imh & x>=1 & x<=imw);
    x = x(tind);     y = y(tind);  b = b(tind);

    % fill polygon region
    if b=='p'
        [x, y, b] = ginput;
        x = round(x);
        y = round(y);        
        tind = y>=1 & y<=imh & x>=1 & x<=imw;
        x = x(tind);     y = y(tind);  b = b(tind);        
        
        lab = str2num(char(b(1)));

        if lab >=0 && lab <= nc
            sp = get_segs_in_poly(x, y, imh, meanx, meany);
            %sp = sp(labels(sp)==0);
            labels(sp) = lab;
        end        

    else
        lab = str2num(char(b));
        try
            if lab >=0 && lab <= nc                
                sp = segimage(y, x);
                labels(sp) = lab;                
            end
        catch
        end
    end

    tmplim = labels(segimage); 
    hue = tmplim / nc;
    lim = hsv2rgb(hue, hue>0, ones(size(hue)))*0.5 + 0.5*grayim;    
    %tmplim = labels(segimage); 
    %tmplim(1:7) = [1:7];
    %lim = im2double(label2rgb(tmplim, 'jet', [1 1 1]))*0.5 + repmat(grayim, [1 1 3])*0.5;
    figure(1), imagesc(lim), axis image      
end

vnames = {'000', '090', '090', '090', '090', '090', 'sky'};
hnames = {'---', '045', '090', '135', 'por', 'sol', '---'};

l2v = [1 2 2 2 2 2 3];
l2h = [0 1 2 3 4 5 0];
ind = labels>0;
imsegs.labels = labels;
imsegs.vert_labels(ind) = l2v(labels(ind));
imsegs.vlabels = repmat({'---'}, [nseg 1]);
ind2 = imsegs.vert_labels>0;
imsegs.vlabels(ind2) = vnames(imsegs.vert_labels(ind2));
imsegs.horz_labels(ind) = l2h(labels(ind));
imsegs.hlabels = repmat({'---'}, [nseg 1]);
ind2 = imsegs.horz_labels>0;
imsegs.hlabels(ind2) = hnames(imsegs.horz_labels(ind2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5



function s = get_segs_in_poly(px, py, height, meanx, meany)       
s = inpolygon(meanx, meany, px, py);
s = find(s);

