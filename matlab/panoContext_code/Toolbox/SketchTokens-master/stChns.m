function chns = stChns( I, opts )
% Compute channels for sketch token detection.
%
% USAGE
%  chns = stChns( I, opts )
%
% INPUTS
%  I          - [h x w x 3] color input image
%  opts       - sketch token model options
%
% OUTPUTS
%  chns       - [h x w x nChannel] output channels
%
% EXAMPLE
%
% See also stDetect, gradientMag, gradientHist
%
% Sketch Token Toolbox     V0.95
% Copyright 2013 Joseph Lim [lim@csail.mit.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License [see bsd.txt]

    % extract gradient magnitude and histogram channels
    if isfield(opts, 'inputColorChannel') && strcmp(opts.inputColorChannel, 'luv')    
        I = rgbConvert(I,'orig');
    else
        I = rgbConvert(I,'luv');
    end
    chns=cell(1,1000);
    k=1;
    chns{k}=I;
    for i = 1:length( opts.sigmas )
        if ( opts.sigmas(i)==0 ),
            I1=I;
        else
            f = fspecial('gaussian',opts.radius,opts.sigmas(i));
            I1 = imfilter(I,f);
        end
        if( opts.nOrients(i)>0 )
            [M,O] = gradientMag( I1, 0, opts.normRad, opts.normConst );
            H = gradientHist( M, O, 1, opts.nOrients(i), 0 );
            k=k+1;
            chns{k}=M;
            k=k+1;
            chns{k}=H;
        else
            M = gradientMag( I1, 0, opts.normRad, opts.normConst );
            k=k+1;
            chns{k}=M;
        end
    end
    chns=cat(3,chns{1:k});
    chns=convTri(chns,opts.chnsSmooth);
end
