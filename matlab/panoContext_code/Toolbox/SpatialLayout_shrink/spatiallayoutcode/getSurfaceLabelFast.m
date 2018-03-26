function [ cimages, Polyg ] = getSurfaceLabelFast(imdir,imagename,workspcdir, layout, vpdata)

moveFolder = './roomModel/SpatialLayout_shrink/';

outimgdir = [workspcdir 'Images/'];
if ~exist(outimgdir,'dir')
    mkdir(outimgdir);
end
workspcdir = [workspcdir 'data/'];
if ~exist(workspcdir,'dir')
    mkdir(workspcdir);
end

% tempdir='./sptiallayouttempworkspace/';
%Compute vps
boxlayout=[];
surface_labels=[];

img=imread([imdir imagename]);
[h w kk]=size(img);
% if isempty(vpdata)
%     [vp p All_lines]=getVP(imdir,imagename,0,workspcdir);
%     VP=vp;
%     if numel(VP)<6
%         aa = zeros(h,w);
%         indd = zeros(h,w);
%         return;
%     end
% 
% 
%     vp=[VP(1) VP(2);VP(3) VP(4);VP(5) VP(6)];
%     [vp P]=ordervp(vp,h,w,p);
%     [vv linemem]=max(P,[],2);
%     vpdata.vp=vp;
%     vpdata.lines=All_lines;
%     vpdata.linemem=linemem;
%     vpdata.dim=[h w];
% %     visvp(vpdata,img)
% end


%Get segmentation and GC surface confidence maps
% sigma=num2str(0.8);
% k1=num2str(100);
% min1=num2str(100);
% inputim=[workspcdir imagename(1:end-4) '.ppm'];
% outputim=[workspcdir imagename(1:end-4) '.pnm'];
% [s w]=system(['./segment' ' ' sigma ' '  k1  ' ' min1 ' ' inputim ' ' outputim ]);



segext='pnm';
nsegments=[5 15 25 35 40 60 80 100];
% fn=['../Imsegs/' imagename(1:end-4) '.' segext];
fn=[moveFolder '/Imsegs/' imagename(1:end-4) '.' segext];
imseg = processSuperpixelImage(fn);

tic
imdata = mcmcComputeImageData(im2double(img), imseg);% made changes here
toc


load(fullfile([moveFolder '/LabelClassifier/'], 'Classifiers_gc.mat'));

spfeatures = mcmcGetAllSuperpixelData(imdir, imseg);
[efeatures, adjlist] = mcmcGetAllEdgeData(spfeatures, imseg(1));


nsegments=[5 15 25 35 40 60 80 100];

pE{1} = test_boosted_dt_mc(eclassifier, efeatures{1});
pE{1} = 1 ./ (1+exp(ecal(1)*pE{1}+ecal(2)));
smaps{1} = msCreateMultipleSegmentations(pE{1}, adjlist{1}, ...
    imseg(1).nseg, nsegments);


for k = 1:numel(nsegments)
    if max(smaps{1}(:, k))>0
        segfeatures{1, k} = mcmcGetSegmentFeatures(imseg, ...
            spfeatures{1}, imdata, smaps{1}(:, k), (1:max(smaps{1}(:, k))));
    end
end


%Get surface label confidences initial from GC
normalize = 1;
pg=zeros(imseg.nseg,7);%7 labels
%get P(L/I)
pg = msTest(imseg, segfeatures, smaps, ...
    labelclassifier, segclassifier,normalize);

filename=fullfile(workspcdir, [imagename(1:end-4) '_lc_gc.mat' ]);
save(filename,'pg');

% rl=[255  255 0  0 0 255 0];
% gl=[255 0 255  255 0 0 0];
% bl=[0  0 255 0 255 255 0];
% rl1=[255 202 132 255 0 255 0 191 0 0];
% gl1=[236 255 112 255 191 0 255 21 0 0];
% bl1=[139 112 255  255 255 0 0 133 255 0];
% cimages = msPg2confidenceImages(imseg,pg);
% [aa indd]=max(cimages{1}(:,:,1:6),[],3);
% 
% figure(102);
% clear mask_color;
% mask_r = rl(indd);
% mask_g = gl(indd);
% mask_b = bl(indd);
% mask_color(:,:,1) = mask_r;
% mask_color(:,:,2) = mask_g;
% mask_color(:,:,3) = mask_b;
% 
% hsvmask=rgb2hsv(mask_color);
% hsvmask(:,:,3)=aa*255;
% %     hsvmask(:,:,2)=aa;
% mask_color=hsv2rgb(hsvmask);
% 
% tempimg = double(img)*0.5 + mask_color*0.5;
% %         tempimg =  mask_color;
% subplot(1,2,1); imshow(uint8(tempimg));


%visualize
[h w kk]=size(img);
if isempty(vpdata)
    [vp p All_lines]=getVP(imdir,imagename,0,workspcdir);
    VP=vp;
    if numel(VP)<6
        aa = zeros(h,w);
        indd = zeros(h,w);
        return;
    end


    vp=[VP(1) VP(2);VP(3) VP(4);VP(5) VP(6)];
    [vp P]=ordervp(vp,h,w,p);
    [vv linemem]=max(P,[],2);
    vpdata.vp=vp;
    vpdata.lines=All_lines;
    vpdata.linemem=linemem;
    vpdata.dim=[h w];
%     visvp(vpdata,img)
end





%Compute intergral images for features
tic
[integData]=getIntegralimages([imagename],vpdata,imseg,500,workspcdir,imdir);
toc



if isempty(layout)
    %Get candidate layouts and their features
    [polyg, Features] = getcandboxlayout( vpdata.vp,vpdata.dim(1),vpdata.dim(2),integData);
else
    [ polyg, Features] = getcandboxlayoutSingle( layout, vpdata.vp,vpdata.dim(1),vpdata.dim(2), integData );
end




%Get initial estimate
Features1=Features;

% load ../LearntClassifiers/pf_i.mat
load([moveFolder '/LearntClassifiers/pf_i.mat']);

lay_score=[];
numL=size(Features,1);
% change features
if(numel(weights)==14)
    tmpf1=sum(Features(:,1:5).*Features(:,11:15),2);
    tmpf2=sum(Features(:,1:5).*Features(:,16:20),2);
    tmpf3=sum(Features(:,6:10).*Features(:,11:15),2);
    tmpf4=sum(Features(:,6:10).*Features(:,16:20),2);
    Features=[Features(:,1:10) tmpf1 tmpf2 tmpf3 tmpf4];
    
end
%evaluate
score=repmat(weights,[numL,1]).*Features;
score=sum(score,2);
[vv ii]=sort(score,'descend');


boxlayout.polyg=polyg;
boxlayout.init=[vv ii];



% load ../LearntClassifiers/pf_il.mat
load([moveFolder '/LearntClassifiers/pf_il.mat']);
Features=Features1;
if(numel(weights)==59)
    tmpf1=sum(Features(:,1:5).*Features(:,11:15),2);
    tmpf2=sum(Features(:,1:5).*Features(:,16:20),2);
    tmpf3=sum(Features(:,6:10).*Features(:,11:15),2);
    tmpf4=sum(Features(:,6:10).*Features(:,16:20),2);
    Features=[Features(:,1:10) tmpf1 tmpf2 tmpf3 tmpf4 Features(:,31:75) ];%Features(:,114:122)];
    
end
score=repmat(weights,[numL,1]).*Features;
score=sum(score,2);
[vv ii]=sort(score,'descend');


boxlayout.reestimated=[vv ii];


lay_scores=[vv ii];
save([workspcdir imagename(1:end-4) '_layres.mat'],'polyg','lay_scores');%,'avg_pg');


% figure(101);
% drawnow;
% for lay=1:25
%     layoutid=ii(lay);
%     Polyg=[];
%     for fie=1:5
%         Polyg{fie}=[];
%         if size(polyg{layoutid,fie})>0
%             Polyg{fie}=polyg{layoutid,fie};
%         end
%     end
%     
%     tempimg=displayout(Polyg,w,h,img);
%     subplot(5,5,lay);imshow(uint8(tempimg),[]);title(num2str(vv(lay)));
% end
% saveas(101,[outimgdir imagename(1:end-4) '_boxlayouts.png']);

Polyg=[];
for fie=1:5
    Polyg{fie}=[];
    if size(polyg{ii(1),fie})>0
        Polyg{fie}=polyg{ii(1),fie};
    end
end
% ShowGTPolyg(img,Polyg,103);
% saveas(103,[outimgdir imagename(1:end-4) '_boxlayout.png']);


%Re-Compute Surface labels (GC+box layout features)
% load(['../LabelClassifier/' 'Classifiers_stage2.mat']);
load([moveFolder '/LabelClassifier/' 'Classifiers_stage2.mat']);
tic
xspfields=[];%save per image
numSup=imseg.nseg;

polyg = Polyg;
numL = size(polyg, 1);

for lay=1:numL %each layout
    xspf=[];
    for  supno=1:numSup
        tempbndy=imdata.tracedbndy{supno}{1};
        if size(tempbndy,1) > 30
            
            YY1=tempbndy(1:10:end,1);XX1=tempbndy(1:10:end,2);
        else
            
            YY1=tempbndy(:,1);XX1=tempbndy(:,2);
        end
        
        for fi=1:5  %each field
            xarea=0;
            if size(polyg{lay,fi},1)>0
                
                XX2=polyg{lay,fi}(:,1);
                YY2=polyg{lay,fi}(:,2);
                
                [in on]=inpolygon(XX1,YY1,[XX2;XX2(1)],[YY2;YY2(1)]);
                
                if numel(find(in==1))==length(in)
                    X0=XX1;Y0=YY1;
                    xarea=polyarea([X0;X0(1)],[Y0;Y0(1)]);
                elseif  numel(find(in==1))==0
                    
                    X0=[];Y0=[];
                    xarea=0;
                else
                    xarea=polyintarea(XX1,YY1,XX2,YY2,0);
                    
                end
                
                
            end
            xspf(supno,fi)=xarea;
        end
        
    end
    xspfields{lay}=xspf;
end
toc



nsp=imseg.nseg; %num of suppixels
smap = [1:nsp];
smaps{1} = smap(:);



clear segfeatures
for k = 1: 1%numel(nsegments)
    features = mcmcGetSegmentFeatures(imseg, ...
        spfeatures{1}, imdata, smaps{1}(:, k), (1:max(smaps{1}(:, k))));
    
    
    
    for lay=1:numL
        tempfeatures=features;
        [fieldfeatures]=segfieldsfeat_sup(imseg,imseg.nseg,smaps{1}(:,k),xspfields{lay});
        tempfeatures(:,95:100)=fieldfeatures(1:size(features,1),:);
        tempfeatures(:,101:106)=pg{1}(:,1:6);
        segfeatures{lay,k}  = tempfeatures;
    end
end




% if numel(vv)>100
%     A = [vv(1) 1;vv(100) 1];
% else
%     A = [vv(1) 1;vv(end) 1];
% end
% b = [0;log(49)];
% params = inv(A)*b;
% 
% %     conf = 1./(1+exp(params(1)*vv(1:5)+params(2)));
% conf = 1./(1+exp(params(1)*vv(1)+params(2)));
% conf=conf./sum(conf);
conf = 1;

avg_pg=zeros(imseg.nseg,7);%7 labels

tic
for j=1:1 %take best 5 scored layouts
%     lay=ii(j);
    lay = j;
    %get P(L/F,I)
    pg = msTest(imseg, segfeatures(lay, :), smaps, ...
        labelclassifier);%, segclassifier,normalize);
    avg_pg=avg_pg+pg{1}.*conf(j);
    
    
end
toc

filename=fullfile(workspcdir,[imagename(1:end-4) '_lc_st2.mat' ]);
save(filename,'avg_pg');

surface_labels.restimated={avg_pg};
surface_labels.init=pg;


% rl=[255  255 0  0 0 255 0];
% gl=[255 0 255  255 0 0 0];
% bl=[0  0 255 0 255 255 0];
% rl1=[255 202 132 255 0 255 0 191 0 0];
% gl1=[236 255 112 255 191 0 255 21 0 0];
% bl1=[139 112 255  255 255 0 0 133 255 0];
cimages = msPg2confidenceImages(imseg,{avg_pg});
% [aa indd]=max(cimages{1}(:,:,1:6),[],3);
% 
% figure(102);
% clear mask_color;
% mask_r = rl(indd);
% mask_g = gl(indd);
% mask_b = bl(indd);
% mask_color(:,:,1) = mask_r;
% mask_color(:,:,2) = mask_g;
% mask_color(:,:,3) = mask_b;
% 
% hsvmask=rgb2hsv(mask_color);
% hsvmask(:,:,3)=aa*255;
% %     hsvmask(:,:,2)=aa;
% mask_color=hsv2rgb(hsvmask);
% 
% tempimg = double(img)*0.5 + mask_color*0.5;
% %         tempimg =  mask_color;
% subplot(1,2,2); imshow(uint8(tempimg));
% % saveas(102,[outimgdir imagename(1:end-4) '_surfacelabels.png']);

end

