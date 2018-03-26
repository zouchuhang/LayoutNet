function model = stTrain( varargin )
% Train SketchTokens classifier.
%
% See stDemo for a full demo that include both traininga and application.
%
% Pre-trained models can be downloaded from:
%  http://people.csail.mit.edu/lim/lzd_cvpr2013/st_data.tgz
%
% Please cite the following paper if you end up using the code:
%  Joseph J. Lim, C. Lawrence Zitnick, and Piotr Dollar. "Sketch Tokens: A
%  Learned Mid-level Representation for Contour and and Object Detection,"
%  CVPR2013.
%
% Note: There is a patent pending on the ideas presented in this work so
% this code should only be used for academic purposes.
%
% USAGE
%  model = stTrain( opts )
%
% INPUTS
%  opts       - parameters (struct or name/value pairs)
%   (1) parameters for model and data:
%   .nClusters  - [150] number of clusters to train with
%   .nTrees     - [25] number of trees in forest to train
%   .radius     - [17] radius of sketch token patches
%   .nPos       - [1000] number of positive patches per cluster
%   .nNeg       - [800] number of negative patches per image
%   .negDist    - [2] distance from closest contour defining a negative
%   .minCount   - [4] minimum number of training examples per node
%   (2) parameters for features:
%   .nCells     - [5] number of self similarity cells
%   .normRad    - [5] normalization radius (see gradientMag)
%   .normConst  - [.01] normalization constant (see gradientMag)
%   .nOrients   - [4 4 0] number of orientations for each channel set
%   .sigmas     - [0 1.5 5] gaussian blur for each channel set
%   .chnsSmooth - [2] radius for channel smoothing (using convTri)
%   .fracFtrs   - [1] fraction of features to use to train each tree
%   (3) other parameters:
%   .seed       - [1] seed for random stream (for reproducibility)
%   .modelDir   - ['models/'] target directory for storing models
%   .modelFnm   - ['model'] model filename
%   .clusterFnm - ['clusters.mat'] file containing cluster info
%   .bsdsDir    - ['BSR/BSDS500/data/'] location of BSDS dataset
%
% OUTPUTS
%  model      - trained sketch token detector w the following fields
%   .trees      - learned forest model struct array (see forestTrain)
%   .opts       - input parameters and constants
%   .clusters   - actual cluster centers used to learn tokens
%
% EXAMPLE
%
% See also stGetPatches, stDetect, forestTrain, chnsCompute, gradientMag
%
% Sketch Token Toolbox     V0.95
% Copyright 2013 Joseph Lim [lim@csail.mit.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License [see bsd.txt]

    % get default parameters
    dfs={'nClusters',150, 'nTrees',25, 'radius',17, 'nPos',1000, 'nNeg',800,...
        'negDist',2, 'minCount',4, 'nCells',5, 'normRad',5, 'normConst',0.01, ...
        'nOrients',[4 4 0], 'sigmas',[0 1.5 5], 'chnsSmooth',2, 'fracFtrs',1, ...
        'seed',1, 'modelDir','models/', 'modelFnm','model', ...
        'clusterFnm','clusters.mat', 'bsdsDir','BSR/BSDS500/data/'};
    opts = getPrmDflt(varargin,dfs,1);

    % if forest exists load it and return
    cd(fileparts(mfilename('fullpath')));
    forestDir = [opts.modelDir '/forest/'];
    forestFn = [forestDir opts.modelFnm];
    if exist([forestFn '.mat'], 'file')
        load([forestFn '.mat']);
        return;
    end

    % compute constants and store in opts
    nTrees=opts.nTrees;
    nCells=opts.nCells;
    
    patchSiz=opts.radius*2+1;
    opts.patchSiz=patchSiz;
    
    nChns = size(stChns(ones(2,2,3),opts),3);
    opts.nChns=nChns;
    
    opts.nChnFtrs = patchSiz*patchSiz*nChns;
    opts.nSimFtrs = (nCells*nCells)*(nCells*nCells-1)/2*nChns;
    opts.nTotFtrs = opts.nChnFtrs + opts.nSimFtrs;
    opts.cellRad = round(patchSiz/nCells/2);
    tmp=opts.cellRad*2+1;
    opts.cellStep = tmp-ceil((nCells*tmp-patchSiz)/(nCells-1)); disp(opts);
    assert( mod(nCells,2)==1 && (nCells-1)*opts.cellStep+tmp <= patchSiz );

    % generate stream for reproducibility of model
    stream=RandStream('mrg32k3a','Seed',opts.seed);

    % train nTrees random trees (can be trained with parfor if enough memory)
    for i=1:nTrees
        stTrainTree( opts, stream, i );
    end

    % accumulate trees and merge into final model
    treeFn = [opts.modelDir '/tree/' opts.modelFnm '_tree'];
    for i=1:nTrees
        t=load([treeFn int2str2(i,3) '.mat'],'tree');
        t=t.tree;
        if (i==1)
            trees=t(ones(1,nTrees));
        else
            trees(i)=t;
        end
    end
    nNodes=0;
    for i=1:nTrees
        nNodes=max(nNodes,size(trees(i).fids,1));
    end
    model.thrs=zeros(nNodes,nTrees,'single');
    Z=zeros(nNodes,nTrees,'uint32');
    model.fids=Z;
    model.child=Z;
    model.count=Z;
    model.depth=Z;
    model.distr=zeros(nNodes,size(trees(1).distr,2),nTrees,'single');
    for i=1:nTrees, tree=trees(i); nNodes1=size(tree.fids,1);
        model.fids(1:nNodes1,i) = tree.fids;
        model.thrs(1:nNodes1,i) = tree.thrs;
        model.child(1:nNodes1,i) = tree.child;
        model.distr(1:nNodes1,:,i) = tree.distr;
        model.count(1:nNodes1,i) = tree.count;
        model.depth(1:nNodes1,i) = tree.depth;
    end
    model.distr = permute(model.distr, [2 1 3]);
    
    clusters=load(opts.clusterFnm);
    clusters=clusters.clusters;
    
    model.opts = opts;
    model.clusters=clusters.clusters;
    if ~exist(forestDir,'dir')
        mkdir(forestDir);
    end
    save([forestFn '.mat'], 'model', '-v7.3');

end

function stTrainTree( opts, stream, treeInd )
% Train a single tree in forest model.

    % location of ground truth
    trnImgDir = [opts.bsdsDir '/images/train/'];
    trnGtDir = [opts.bsdsDir '/groundTruth/train/'];
    imgIds=dir([trnImgDir '*.jpg']);
    imgIds={imgIds.name};
    nImgs=length(imgIds);
    for i=1:nImgs,
        imgIds{i}=imgIds{i}(1:end-4);
    end

    % extract commonly used options
    radius=opts.radius;
    patchSiz=opts.patchSiz;
    nChns=opts.nChns;
    nTotFtrs=opts.nTotFtrs;
    nClusters=opts.nClusters;
    nPos=opts.nPos;
    nNeg=opts.nNeg;

    % finalize setup
    treeDir = [opts.modelDir '/tree/'];
    treeFn = [treeDir opts.modelFnm '_tree'];
    if exist([treeFn int2str2(treeInd,3) '.mat'],'file')
        return;
    end
    fprintf('\n-------------------------------------------\n');
    fprintf('Training tree %d of %d\n',treeInd,opts.nTrees);
    tStart=clock;

    % set global stream to stream with given substream (will undo at end)
    streamOrig = RandStream.getGlobalStream();
    set(stream,'Substream',treeInd);
    RandStream.setGlobalStream( stream );

    % sample nPos positive patch locations per cluster
    clstr=load(opts.clusterFnm);
    clstr=clstr.clusters;
    for i = 1:nClusters
        if i==1
            centers=[];
        end
        ids = find(clstr.clusterId == i);
        ids = ids(randperm(length(ids),min(nPos,length(ids))));
        centers = [centers; [clstr.x(ids),clstr.y(ids),clstr.imId(ids),...
            clstr.clusterId(ids),clstr.gtId(ids)]]; %#ok<AGROW>
    end

    % collect positive and negative patches and compute features
    fids=sort(randperm(nTotFtrs,round(nTotFtrs*opts.fracFtrs)));
    k = size(centers,1)+nNeg*nImgs;
    ftrs = zeros(k,length(fids),'single');
    labels = zeros(k,1); k = 0;
    tid = ticStatus('Collecting data',1,1);
    for i = 1:nImgs
        % get image and compute channels
        gt=load([trnGtDir imgIds{i} '.mat']);
        gt=gt.groundTruth;
        
        I = imread([trnImgDir imgIds{i} '.jpg']);
        I = imPad(I,radius,'symmetric');
        chns = stChns(I,opts);
        
        % sample positive patch locations
        centers1=centers(centers(:,3)==i,:);
        lbls1=centers1(:,4);
        xy1=single(centers1(:,[1 2]));
        
        % sample negative patch locations
        M=false(size(I,1)-2*radius,size(I,2)-2*radius);
        nGt=length(gt);
        for j=1:nGt
            M1=gt{j}.Boundaries;
            if ~isempty(M1)
                M=M | M1;
            end
        end
        M(bwdist(M)<opts.negDist)=1;
        M=~M;
        M([1:radius end-radius:end],:)=0;
        M(:,[1:radius end-radius:end])=0;
        [y,x]=find(M);
        k1=min(length(y),nNeg);
        
        rp=randperm(length(y),k1);
        y=y(rp);
        x=x(rp);
        xy0=[x y];
        lbls0=ones(k1,1)*(nClusters+1);
        
        % crop patches
        xy=[xy1; xy0];
        lbls=[lbls1; lbls0];
        k1=length(lbls);
        ps=zeros(patchSiz,patchSiz,nChns,k1,'single');
        p=patchSiz-1;
        for j=1:k1
            ps(:,:,:,j)=chns(xy(j,2):xy(j,2)+p,xy(j,1):xy(j,1)+p,:);
        end
        if(0), montage2(squeeze(ps(:,:,1,:))); drawnow; end
        
        % compute features and store
        ftrs1=[reshape(ps,[],k1)' stComputeSimFtrs(ps,opts)];
        ftrs(k+1:k+k1,:) = ftrs1(:,fids);
        labels(k+1:k+k1) = lbls;
        k=k+k1;
        
        tocStatus(tid,i/nImgs);
    end
    if k<size(ftrs,1)
        ftrs=ftrs(1:k,:);
        labels=labels(1:k);
    end

    % train sketch token classifier (random decision tree)
    tree=forestTrain(ftrs,labels,'maxDepth',999);
    tree.fids(tree.child>0) = fids(tree.fids(tree.child>0)+1)-1;
    tree=pruneTree(tree,opts.minCount); %#ok<NASGU>
    if ~exist(treeDir,'dir')
        mkdir(treeDir);
    end
    save([treeFn int2str2(treeInd,3) '.mat'],'tree');
    e=etime(clock,tStart);
    fprintf('Training of tree %d complete (time=%.1fs).\n',treeInd,e);
    RandStream.setGlobalStream( streamOrig );

end

function tree = pruneTree( tree, minCount )
% Prune all nodes whose count is less than minCount.

    % mark all internal nodes if either child has count<=minCount
    mark = [0; tree.count<=minCount];
    mark = mark(tree.child+1) | mark(tree.child+2);
    
    % list of nodes to be discarded / kept
    disc=tree.child(mark);
    disc=[disc; disc+1];
    n=length(tree.fids);
    keep=1:n;
    keep(disc)=[];
    
    % prune tree
    tree.fids=tree.fids(keep);
    tree.thrs=tree.thrs(keep);
    tree.child=tree.child(keep);
    tree.distr=tree.distr(keep,:);
    tree.count=tree.count(keep);
    tree.depth=tree.depth(keep);
    assert(all(tree.count>minCount))
    
    % re-index children
    route=zeros(1,n);
    route(keep)=1:length(keep);
    tree.child(tree.child>0) = route(tree.child(tree.child>0));
end

function ftrs = stComputeSimFtrs( chns, opts )
% Compute self-similarity features.
    n=opts.nCells;
    if(n==0),
        ftrs=[];
        return;
    end
    nSimFtrs=opts.nSimFtrs;
    nChns=opts.nChns;
    m=size(chns,4);
    
    inds = ((1:n)-(n+1)/2)*opts.cellStep+opts.radius+1;
    chns=reshape(chns,opts.patchSiz,opts.patchSiz,nChns*m);
    chns=convBox(chns,opts.cellRad);
    chns=reshape(chns(inds,inds,:,:),n*n,nChns,m);
    ftrs=zeros(nSimFtrs/nChns,nChns,m,'single');
    k=0;
    for i=1:n*n-1
        k1=n*n-i;
        ftrs(k+1:k+k1,:,:)=chns(1:end-i,:,:)-chns(i+1:end,:,:);
        k=k+k1;
    end
    ftrs = reshape(ftrs,nSimFtrs,m)';
    % % For m=1, the above should be identical to the following:
    % [cids1,cids2]=computeCids(size(chns),opts); % see stDetect.m
    % chns=convBox(chns,opts.cellRad); k=opts.nChnFtrs;
    % cids1=cids1(k+1:end)-k+1; cids2=cids2(k+1:end)-k+1;
    % ftrs=chns(cids1)-chns(cids2);
end
