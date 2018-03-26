function hBlobs = RecreateBlobHierarchy(blobs, hierarchy)
% [blobs hierarchy] = RecreateBlobHierarchy(blobs, hierarchy)
% 
% Recreates the hierarchical grouping using the starting blobs and the 
% resulting hierarchy. This allows one to save the grouping using
% relatively small disk space while still being able to fastly recreate the
% complete grouping.
%
% blobs:            Input cell array with blobs
% hierarchy:        Hierarchy of the blobs as created by
%                   HierarchicalGrouping.m
%
% hBlobs:           All segments of the hierarchical grouping.
%
%     Jasper Uijlings - 2013

hBlobs = cell(length(hierarchy) + 1,1);

hBlobs(1:length(blobs)) = blobs;

for i=length(blobs)+1:length(hBlobs)
    n = find(hierarchy == i);
    
    if length(n) ~= 2
        error('One can not merge more than 2 blobs!');
    end
    
    hBlobs{i} = MergeBlobs(hBlobs{n(1)}, hBlobs{n(2)});
end