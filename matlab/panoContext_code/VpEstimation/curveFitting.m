function [ outputXYZ, outputNM ] = curveFitting( inputXYZ, weight )
%CURVEFITTING Summary of this function goes here
%   Detailed explanation goes here
len = sqrt(sum(inputXYZ.^2,2));
inputXYZ = inputXYZ ./ repmat(len,[1 3]);
weightXYZ = inputXYZ .* repmat(weight,[1 3]);
XX = sum(weightXYZ(:,1).^2);
YY = sum(weightXYZ(:,2).^2);
ZZ = sum(weightXYZ(:,3).^2);
XY = sum(weightXYZ(:,1).*weightXYZ(:,2));
YZ = sum(weightXYZ(:,2).*weightXYZ(:,3));
ZX = sum(weightXYZ(:,3).*weightXYZ(:,1));

A = [XX XY ZX; ...
     XY YY YZ; ...
     ZX YZ ZZ];
[~,~,V] =  svd(A);
outputNM = V(:,end)';
outputNM = outputNM./norm(outputNM);
outputXYZ = [];

end

