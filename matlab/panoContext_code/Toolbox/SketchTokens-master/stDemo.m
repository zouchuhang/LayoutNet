%% Sketch Token demo and usage.
%%
%% Please cite the following paper if you end up using the code:
%%  Joseph J. Lim, C. Lawrence Zitnick, and Piotr Dollar. "Sketch Tokens: A
%%  Learned Mid-level Representation for Contour and and Object Detection,"
%%  CVPR2013.
%%
%% Note: There is a patent pending on the ideas presented in this work so
%% this code should only be used for academic purposes.
%%
%% Sketch Token Toolbox     V0.95
%% Copyright 2013 Joseph Lim [lim@csail.mit.edu]
%% Please email me if you find bugs, or have suggestions or questions!
%% Licensed under the Simplified BSD License [see bsd.txt]



%% setup (follow instructions, only need to do once)
cd(fileparts(mfilename('fullpath')))
if( 0 )
    % (1) Download Berkeley Segmentation Data Set and Benchmarks 500
    % http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/
    % BSR/ should contain BSDS500, bench, and documentation
    addpath('BSR/bench/benchmarks');
    % (2) Download and compile Piotr's toolbox (directions on website)
    % http://vision.ucsd.edu/~pdollar/toolbox/doc/index.html
    % (3) Compile code (removed OPTIMFLAGS to compile single core variant)
    mex stDetectMex.cpp 'OPTIMFLAGS="$OPTIMFLAGS' '/openmp"'
    % (4) Add current directory to path
    addpath(pwd);
end

%% train or load model (see stTrain.m)
if( 0 )
    if( 1 )
        % can be trained in ~20m and requires ~3GB ram
        % BSDS500 performance: ODS=0.721 OIS=0.739 AP=0.768
        opts=struct('nPos',100,'nNeg',80,'modelFnm','modelSmall','nTrees',20);
    else
        % takes ~2.5 hours to train and requires ~27GB ram
        % BSDS500 performance: ODS=0.727 OIS=0.746 AP=0.780
        opts=struct('nPos',1000,'nNeg',800,'modelFnm','modelFull');
    end
    tic, model=stTrain(opts); toc
else
    % Pre-trained models can be downloaded from:
    %  http://people.csail.mit.edu/lim/lzd_cvpr2013/st_data.tgz
    load('models/forest/modelSmall.mat');
end

%% evaluate sketch token model on BSDS500 (see stEvalBsds.m)
if (0),
    [ODS,OIS,AP]=stEvalBsds( model );
end

%% detect sketch tokens and extract edges (see stDetect.m and stToEdges.m)
I = imread('peppers.png');
tic; st = stDetect( I, model ); toc;
tic, E = stToEdges( st, 1 ); toc

%% visualize edge detection results
figure(1);
im(I);
figure(2);
im(E);
colormap jet;
