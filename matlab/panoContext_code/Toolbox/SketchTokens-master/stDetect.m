function S = stDetect( I, model, stride, rescale_back )
% Detect sketch tokens in image.
%
% USAGE
%  S = stDetect( I, model, [stride] )
%
% INPUTS
%  I            - [h x w x 3] color input image
%  model        - sketch token model trained with stTrain
%  stride       - [2] stride at which to compute sketch tokens
%  rescale_back - [true] rescale after running stride
%
% OUTPUTS
%  S          - [h x w x (nTokens+1)] sketch token probability maps
%
% EXAMPLE
%
% See also stTrain, stChns
%
% Sketch Token Toolbox     V0.95
% Copyright 2013 Joseph Lim [lim@csail.mit.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License [see bsd.txt]

    if nargin<3
        stride=2;
    end
    if nargin<4
        rescale_back = true;
    end

    % compute features
    sizeOrig=size(I);
    opts=model.opts;
    %opts.inputColorChannel = 'luv';
    opts.inputColorChannel = 'rgb';
    I = imPad(I,opts.radius,'symmetric');
    chns = stChns( I, opts );
    [cids1,cids2] = computeCids(size(chns),opts);
    chnsSs = convBox(chns,opts.cellRad);

    % run forest on image
    S = stDetectMex( chns, chnsSs, model.thrs, model.fids, model.child, ...
      model.distr, cids1, cids2, stride, opts.radius, opts.nChnFtrs );

    % finalize sketch token probability maps
    S = permute(S,[2 3 1]) * (1/opts.nTrees);
    if ~rescale_back
        %keyboard;
    else    
        S = imResample( S, stride );
        cr=size(S); cr=cr(1:2)-sizeOrig(1:2);
        if any(cr)
            S=S(1:end-cr(1),1:end-cr(2),:);
        end
    end

end

function [cids1,cids2] = computeCids( siz, opts )
    % construct cids lookup for standard features
    radius=opts.radius;
    s=opts.patchSiz;
    nChns=opts.nChns;
    
    ht=siz(1);
    wd=siz(2);
    assert(siz(3)==nChns);
    
    nChnFtrs=s*s*nChns;
    fids=uint32(0:nChnFtrs-1);
    rs=mod(fids,s);
    fids=(fids-rs)/s;
    cs=mod(fids,s);
    ch=(fids-cs)/s;
    cids = rs + cs*ht + ch*ht*wd;
    
    % construct cids1/cids2 lookup for self-similarity features
    n=opts.nCells;
    m=opts.cellStep;
    nCellTotal=(n*n)*(n*n-1)/2;
    
    assert(mod(n,2)==1); n1=(n-1)/2;
    nSimFtrs=nCellTotal*nChns;
    fids=uint32(0:nSimFtrs-1);
    ind=mod(fids,nCellTotal);
    ch=(fids-ind)/nCellTotal;
    
    k=0;
    for i=1:n*n-1,
        k1=n*n-i;
        ind1(k+1:k+k1)=(0:k1-1);
        k=k+k1;
    end
    k=0;
    for i=1:n*n-1,
        k1=n*n-i;
        ind2(k+1:k+k1)=(0:k1-1)+i;
        k=k+k1;
    end
    
    ind1=ind1(ind+1);
    rs1=mod(ind1,n);
    cs1=(ind1-rs1)/n;
    ind2=ind2(ind+1);
    rs2=mod(ind2,n);
    cs2=(ind2-rs2)/n;
    
    rs1=uint32((rs1-n1)*m+radius);
    cs1=uint32((cs1-n1)*m+radius);
    rs2=uint32((rs2-n1)*m+radius);
    cs2=uint32((cs2-n1)*m+radius);
    
    cids1 = rs1 + cs1*ht + ch*ht*wd;
    cids2 = rs2 + cs2*ht + ch*ht*wd;
    
    % combine cids for standard and self-similarity features
    cids1=[cids cids1];
    cids2=[zeros(1,nChnFtrs) cids2];
end
