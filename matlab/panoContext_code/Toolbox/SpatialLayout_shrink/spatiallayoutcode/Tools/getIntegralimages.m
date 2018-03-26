function [integData]=getIntegralimages(imagename,vpdata,imseg,quantsiz,labelconfdir,imdir)
% getIntegralimages Given an image, and its features compute the integral images of
% rectified features. This code gives integral images for line membership
% features and surface label features.
%
%
% INPUT:
%   imagename - original image name
%    vpdata- detected line segments  and vanishing points
%    imseg, labelconfdir - surface label features
%    quantsiz- size of the rectified image
% OUTPUT:
%    Integral image for each feature rectified in direction of each pair of
%    vanishing point 12, 23 and 31.

%For more details check [2] Varsha Hedau, Derek Hoiem, David Forsyth, Thinking Inside the Box:
%    Using Appearance Models and Context Based on Room Geometry, ECCV 2010.


LineNum=4;
LabelNum=7;


angle=180;
nbins=6;

%surface labels
filename=fullfile(labelconfdir,[imagename(1:end-4) '_lc_gc.mat' ]);
if exist(filename)
    load(filename);
    cimages = msPg2confidenceImages(imseg,pg);
end


h=vpdata.dim(1);
w=vpdata.dim(2);


% line membership features
for l=1:LineNum
    ind=find(vpdata.linemem==l);
    Lines=vpdata.lines(ind,:);
    x1=Lines(:,1);x2=Lines(:,2);y1=Lines(:,3);y2=Lines(:,4);
    tx1=repmat(x1,[1,11]);tx2=repmat(x2,[1,11]);
    ty1=repmat(y1,[1,11]);ty2=repmat(y2,[1,11]);
    tt2=[0:0.1:1];tt1=[1:-0.1:0];
    tt1=repmat(tt1,[size(Lines,1),1]);tt2=repmat(tt2,[size(Lines,1),1]);
    xx=tx2.*tt2+tx1.*tt1;
    yy=ty2.*tt2+ty1.*tt1;
    xx=xx(:);yy=yy(:);
    xx=min([w*ones(size(xx,1)) round(xx)],[],2);
    xx=max([ones(size(xx,1)) round(xx)],[],2);
    yy=min([h*ones(size(yy,1)) round(yy)],[],2);
    yy=max([ones(size(yy,1)) round(yy)],[],2);
    
    tlen=Lines(:,7);
    tlen=repmat(tlen,[1 11]);
    tlen=tlen(:);
    
    feat_img{l}=zeros(h,w);
    indd=sub2ind(size(feat_img{l}),yy,xx);
    feat_img{l}(indd)=tlen/11;
    
end

%Weighted line membership features
objconf=double(cimages{1}(:,:,6));%weighted with object conf
for l=1:LineNum
    feat_img{LineNum+l}=(1-objconf).*feat_img{l};
end

%surface label features
for l=1:LabelNum
    feat_img{2*LineNum+l}=double(cimages{1}(:,:,l));
end



integData=struct('intI12',' ','intI23',' ','intI31',' ',...
    'intnum12',' ', 'intnum23',' ', 'intnum31',' ');


for feat=1: length(feat_img)
    
    intI12{feat}=[];
    intI23{feat}=[];
    intI31{feat}=[];
    
    
    vno1=1;vno2=2;
    [feat_avg, feat_sum, num_12, theta_ind12]= txfmImg(feat_img{feat},vpdata.vp([vno1 vno2],:),0,quantsiz,1);
    intI12{feat} = cumsum(cumsum(double(feat_sum)),2);
    intI12{feat}=[zeros(1,size(intI12{feat},2)) ; intI12{feat}];
    intI12{feat}=[zeros(size(intI12{feat},1),1)  intI12{feat}];
    
    vno1=2;vno2=3;
    [feat_avg, feat_sum, num_23, theta_ind23]= txfmImg(feat_img{feat},vpdata.vp([vno1 vno2],:),0,quantsiz,1);
    intI23{feat} = cumsum(cumsum(double(feat_sum)),2);
    intI23{feat}=[zeros(1,size(intI23{feat},2)) ; intI23{feat}];
    intI23{feat}=[zeros(size(intI23{feat},1),1)  intI23{feat}];
    
    vno1=3;vno2=1;
    [feat_avg, feat_sum, num_31, theta_ind31]= txfmImg(feat_img{feat},vpdata.vp([vno1 vno2],:),0,quantsiz,1);
    intI31{feat} = cumsum(cumsum(double(feat_sum)),2);
    intI31{feat}=[zeros(1,size(intI31{feat},2)) ; intI31{feat}];
    intI31{feat}=[zeros(size(intI31{feat},1),1)  intI31{feat}];
    
    %check for nan inf
    if numel(find(isnan(feat_sum))) > 0 | numel(find(isinf(feat_sum))) > 0
        disp('Feats are NaN');
        keyboard;
    end
    
end



%% save integData
intnum12=cumsum(cumsum(double(num_12)),2);
intnum12=[zeros(1,size(intnum12,2)) ; intnum12];
intnum12=[zeros(size(intnum12,1),1)  intnum12];

intnum23=cumsum(cumsum(double(num_23)),2);
intnum23=[zeros(1,size(intnum23,2)) ; intnum23];
intnum23=[zeros(size(intnum23,1),1)  intnum23];

intnum31=cumsum(cumsum(double(num_31)),2);
intnum31=[zeros(1,size(intnum31,2)) ; intnum31];
intnum31=[zeros(size(intnum31,1),1)  intnum31];


integData.intI12=intI12;
integData.intI23=intI23;
integData.intI31=intI31;
integData.intnum12=intnum12;
integData.intnum23=intnum23;
integData.intnum31=intnum31;

integData.theta_ind12=theta_ind12;
integData.theta_ind23=theta_ind23;
integData.theta_ind31=theta_ind31;
%%


return;