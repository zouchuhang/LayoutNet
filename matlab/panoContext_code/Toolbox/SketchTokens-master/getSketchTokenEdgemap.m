function [ E ] = getSketchTokenEdgemap( img )
%GETSKETCHTOKENEDGEMAP Summary of this function goes here
%   Detailed explanation goes here
load('./Toolbox/SketchTokens-master/models/forest/modelSmall.mat');
st = stDetect( img, model );
E = stToEdges( st, 1 );

end

