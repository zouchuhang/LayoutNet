function [similarity indSim] = SSSimSize(a, b, blobStruct)
% function similarity = SSSimSize(a, b, blobStruct)
%
% Calculate size similarity

similarity = (blobStruct.imSize - blobStruct.size(a) - blobStruct.size(b)) ...
           ./ blobStruct.imSize;

indSim = similarity;