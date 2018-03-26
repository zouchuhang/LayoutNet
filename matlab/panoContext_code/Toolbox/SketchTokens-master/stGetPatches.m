function clusters = stGetPatches( radius, nPatches, bsdsDir )
% Sample ground truth edge sketch patches.
%
% Calling stGetPatches() is *optional* as the code package comes with
% pre-computed clusters.mat. stGetPatches() only needs to be called
% if you want to generate alternate clusters to the ones provided.
%
% stGetPatches() is the first step in generating sketch tokens classes: it
% is used to sample patches of ground truth edge maps from the training set
% of the Berkeley Segmentation Dataset. After sampling, the extracted
% patches must be clustered to generate the sketch token classes.
% Clustering code is *NOT* provided. It should be easy to implement, see
% the paper for more details. After clustering, the additional fields
% "clusterId" and "clusters" should be initialized in the clusters struct.
% "clusterId" should indicate cluster membership of each extracted patch
% (integer between 1 and K) and "cluters" should be the mean of the patches
% belonging to the given cluster. After these two fields are added to the
% clusters struct, the sketch token model is ready to be trained via
% stTrain() (see parameter "clusterFnm" in stTrain.m).
%
% USAGE
%  clusters = stGetPatches( [radius], [nPatches], [bsdsDir] )
%
% INPUTS
%  radius     - [15] radius of sketch token patches
%  nPatches   - [inf] maximum number of patches to sample
%  bsdsDir    - ['BSR/BSDS500/data/'] location of BSDS dataset
%
% OUTPUTS
%  clusters   - extracted ground truth info w the following fields
%   .x          - [Nx1] x-coordinate each sampled patch
%   .y          - [Nx1] y-coordinate each sampled patch
%   .gtId       - [Nx1] integer ground turth labeler of each sampled patch
%   .imId       - [Nx1] integer image id of each sampled patch
%   .patches    - [PxPxN] binary images of sampled patches, P=2*radius+1
%   .clusterId  - [Nx1] cluster membership in [1,K] (NOT COMPUTED)
%   .clusters   - [PxPxK] cluster images (NOT COMPUTED)
%
% EXAMPLE
%
% See also stTrain
%
% Sketch Token Toolbox     V0.95
% Copyright 2013 Joseph Lim [lim@csail.mit.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License [see bsd.txt]

    if( nargin<1 ),
        radius=15;
    end
    if( nargin<2 ),
        nPatches=inf;
    end
    if( nargin<3 ),
        bsdsDir='BSR/BSDS500/data/';
    end

    % location of ground truth
    trnImgDir = [bsdsDir '/images/train/'];
    trnGtDir = [bsdsDir '/groundTruth/train/'];
    imgIds=dir([trnImgDir '*.jpg']);
    imgIds={imgIds.name};
    nImgs=length(imgIds);
    for i=1:nImgs,
        imgIds{i}=imgIds{i}(1:end-4);
    end

    % loop over ground truth and collect samples
    clusters=struct('x',[],'y',[],'gtId',[],'imId',[],'patches',[],...
      'clusterId',[],'clusters',[]);
    clusters.patches = false(radius*2+1,radius*2+1,9000*5*nImgs);
    tid = ticStatus('data collection');
    cnt=0;
    for i = 1:nImgs
        gt=load([trnGtDir imgIds{i} '.mat']);
        gt=gt.groundTruth;
        for j=1:length(gt)
            if(isempty(gt{j}.Boundaries)),
                continue;
            end
            M0 = gt{j}.Boundaries;
            M=M0;
            M([1:radius end-radius+1:end],:)=0;
            M(:,[1:radius end-radius+1:end])=0;
            [y,x]=find(M);
            cnt1=length(y);
            clusters.y = [clusters.y; int32(y)];
            clusters.x = [clusters.x; int32(x)];
            clusters.gtId = [clusters.gtId; ones(cnt1,1,'int32')*j];
            clusters.imId = [clusters.imId; ones(cnt1,1,'int32')*i];
            for k=1:cnt1,
                clusters.patches(:,:,cnt+k) = M0(y(k)-radius:y(k)+radius,x(k)-radius:x(k)+radius);
            end
            cnt = cnt + cnt1;
        end
        tocStatus(tid, i/nImgs);
    end
    clusters.patches = clusters.patches(:,:,1:cnt);

    % optionally sample patches
    if( nPatches<cnt )
        stream = RandStream('mrg32k3a','Seed',1);
        kp = sort(randperm(stream,cnt,nPatches));
        clusters.x=clusters.x(kp);
        clusters.gtId=clusters.gtId(kp);
        clusters.y=clusters.y(kp);
        clusters.imId=clusters.imId(kp);
        clusters.patches=clusters.patches(:,:,kp);
    end

end
