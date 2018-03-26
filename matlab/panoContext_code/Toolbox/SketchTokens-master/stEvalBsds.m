function [ODS,OIS,AP] = stEvalBsds( model, varargin )
% Evaluate sketch token edge detector on BSDS500.
%
% USAGE
%  [ODS,OIS,AP] = stEvalBsds( model, parameters )
%
% INPUTS
%  model      - sketch token model or model filename
%  parameters - parameters (struct or name/value pairs)
%   .dataType   - ['test'] should be either 'test' or 'val'
%   .nThresh    - [99] number of thresholds for evaluation
%   .cleanup    - [0] if true delete temporary files
%   .show       - [0] figure for displaying results (or 0)
%   .modelDir   - [] directory for storing models
%   .bsdsDir    - [] directory of BSDS dataset
%   .name       - [''] name to append to evaluation
%   .pDistr     - [{'type','parfor'}] parameters for fevalDistr
%
% OUTPUTS
%  ODS        - standard error measure on BSDS500
%  OIS        - standard error measure on BSDS500
%  AP         - standard error measure on BSDS500
%
% EXAMPLE
%
% See also stDetect, stTrain
%
% Sketch Token Toolbox     V0.95
% Copyright 2013 Joseph Lim [lim@csail.mit.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License [see bsd.txt]

    % get default parameters
    dfs={'dataType','test', 'nThresh',99, 'cleanup',0, 'show',0, ...
      'modelDir',[], 'bsdsDir',[], 'name','', 'pDistr',{{'type','parfor'}} };
    p=getPrmDflt(varargin,dfs,1);
    if( ischar(model) ),
        model=load(model);
        model=model.model;
    end
    if( isempty(p.modelDir )),
        p.modelDir=model.opts.modelDir;
    end
    if( isempty(p.bsdsDir )),
        p.bsdsDir=model.opts.bsdsDir;
    end
    p.modelDir = [p.modelDir '/' p.dataType '/'];

    % eval on either validation set or test set
    imgDir = [p.bsdsDir '/images/' p.dataType '/'];
    gtDir = [p.bsdsDir '/groundTruth/' p.dataType '/'];
    evalDir = [p.modelDir model.opts.modelFnm p.name '-eval/'];
    resDir = [p.modelDir model.opts.modelFnm p.name '/'];
    assert(exist(imgDir,'dir')==7);
    assert(exist(gtDir,'dir')==7);

    % if evaluation exists collect results and display
    if(exist([evalDir '/eval_bdry.txt'],'file'))
        [ODS,~,~,~,OIS,~,~,AP]=collect_eval_bdry(evalDir);
        fprintf('ODS=%.3f OIS=%.3f AP=%.3f\n',ODS,OIS,AP);
        if( p.show ),
            plot_eval(evalDir,'r');
        end; 
        return;
    end

    % get image ids
    ids=dir([imgDir '*.jpg']);
    ids={ids.name};
    n=length(ids);
    for i=1:n,
        ids{i}=ids{i}(1:end-4);
    end

    % extract sketch tokens and save results
    if(~exist(resDir,'dir')),
        mkdir(resDir);
    end;
    do=false(1,n);
    
    for i=1:n,
        do(i)=~exist([resDir ids{i} '.png'],'file');
    end
    do=find(do);
    m=length(do);
    parfor i=1:m,
        id=ids{do(i)}; %#ok<PFBNS>
        
        I = imread([imgDir id '.jpg']);
        st = max(min(1,stDetect(I,model)),0);
        E = stToEdges(st,1); 
        
        imwrite(uint8(E*255),[resDir id '.png']);
    end

    % perform evaluation on each image (Linux only, slow)
    if(ispc),
        error('Evaluation code runs on Linux ONLY.');
    end
    do=false(1,n);
    jobs=cell(1,n);
    for i=1:n,
        do(i)=~exist([evalDir ids{i} '_ev1.txt'],'file');
    end
    for i=1:n,
        id=ids{i};
        jobs{i}={[resDir id '.png'],...
        [gtDir id '.mat'],[evalDir id '_ev1.txt'],p.nThresh};
    end
    if(~exist(evalDir,'dir')),
        mkdir(evalDir);
    end
    fevalDistr('evaluation_bdry_image',jobs(do),p.pDistr{:});

    % collect results and display
    [ODS,~,~,~,OIS,~,~,AP]=collect_eval_bdry(evalDir);
    fprintf('ODS=%.3f OIS=%.3f AP=%.3f\n',ODS,OIS,AP);
    if( p.show ),
        plot_eval(evalDir,'r');
    end
    if( p.cleanup ),
        delete([evalDir '/*_ev1.txt']);
        delete([evalDir '/eval_bdry_img.txt']);
        delete([evalDir '/eval_bdry_thr.txt']);
        delete([resDir '/*.png']); rmdir(resDir);
    end

end
