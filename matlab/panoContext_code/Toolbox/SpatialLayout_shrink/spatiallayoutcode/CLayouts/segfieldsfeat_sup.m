function [fieldfeatures]=segfieldsfeat(imsegs,nsegments,smap,xspfield,spfind,spdata,features)
sinds=[1:nsegments];
xspfield(isnan(xspfield))=0;

if exist('spfind','var')
    spfind1{1}=spfind{1};%floor
    spfind1{2}=[spfind{2};spfind{3};spfind{4};spfind{5}];% walls+ceiling
    clear spfind;
    spfind=spfind1;
end

%get model for fields
if exist('spfind','var')
    bgmean=[];
    bgmedian=[];
    for f =1:numel(spfind)
        if numel(spfind1{f})==0
            bgmean(f,1:21)=0;
            bgmedian(f,1:21)=0;
            continue;
        end
        sparea = imsegs.npixels(spfind{f});
        npix = sum(sparea);
        spnorm = sparea / npix;
        
        bgmean(f,1:3) = sum( spdata(spfind{f}, 1:3).*repmat(spnorm, [1 3]), 1);
        bgmean(f,4:6) = rgb2hsv(bgmean(f, 1:3));
        bgmean(f,7:21)=sum(spdata(spfind{f}, 15:29) .* repmat(spnorm, [1 15]), 1);
        
        featinds = [1:3];
        for i=1:length(featinds)
            [vv,ii] = sort(spdata(spfind{f},featinds(i)));
            temp = cumsum(spnorm(ii));
            tempfeat1(i) = vv(min(find(temp>0.5)));
        end
        bgmedian(f,1:3) = tempfeat1(:)';
        bgmedian(f,4:6) = rgb2hsv(bgmedian(f, 1:3));
        
        featinds = [15:29];
        for i=1:length(featinds)
            [vv,ii] = sort(spdata(spfind{f},featinds(i)));
            temp = cumsum(spnorm(ii));
            tempfeat2(i) = vv(min(find(temp>0.5)));
        end
        bgmedian(f,7:21) = tempfeat2(:)';
        
    end
end




for k = 1:nsegments
    
    spind = find(smap==sinds(k));
    sparea = imsegs.npixels(spind);
    npix = sum(sparea);
    
    oarea=sum(xspfield(spind,:),1);
    tfeatures = oarea(:)';
    tfeatures=tfeatures/(npix+(npix==0));
    tempfieldfeatures=tfeatures;
    tempfieldfeatures(tfeatures==0)=eps;
    
    entfieldfeat=-1*sum(tfeatures.*log2(tempfieldfeatures));
    tfeatures=[tfeatures entfieldfeat ];
    fieldfeatures(k,:)=tfeatures;
    
    nf=numel(fieldfeatures);
    if exist('spfind','var')
        
        for tempf=1:size(bgmedian,1)
            fieldfeatures(k,nf+(1:5)) = [abs(bgmedian(tempf,1:4)-features(k,1:4)) ...
                norm((features(k,15:29)- bgmedian(tempf,7:21)),2)] ;%r g b h and mean tex res diff
            nf=nf+5;
        end
    end
    
end



end