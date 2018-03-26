function [ meandists ] = pairwiseHyperTest( seedrect, testrects, pairwise, pairtype, seedtype, testtype, data )
%PAIRWISEHYPERTEST Summary of this function goes here
%   Detailed explanation goes here

if size(testrects,1)==0
    meandists = zeros(0,1);
    return;
end

K = min(5, size(pairwise,1)-1);
if K<3
    meandists = zeros(size(testrects,1),1);
    return;
end

if ismember(seedtype, [1 3 4 6 8]) && pairtype==2 % seed type is rectangle, cannot use normalized
    meandists = zeros(size(testrects,1),1);
    return;
end

if ismember(seedtype, [1 3 4 6 8]) && seedtype==testtype && any((seedrect(4:6)-seedrect(1:3))<0.01) % seed type is the same as test, and both are rectangles
    seedangle = data.seedangle;
    testangle = data.testangle;
    
    intest = testangle==seedangle;
    data.normsize = seedrect(4:6) - seedrect(1:3);
    seedpoint = rect2point(seedrect, pairtype, data);
    testpoint = rect2point(testrects(intest,:), pairtype, data);
    querypoint = testpoint - repmat(seedpoint, size(testpoint,1), 1);
    
    if any(intest)
        [~, dists] = annsearch(pairwise', querypoint', K);
        intestdist = mean(dists, 1);

        meandists = zeros(size(testrects,1),1);
        meandists(intest) = intestdist;
    else
        meandists = zeros(size(testrects,1),1);
    end
    return;
end

data.normsize = seedrect(4:6) - seedrect(1:3);
seedpoint = rect2point(seedrect, pairtype, data);
testpoint = rect2point(testrects, pairtype, data);
querypoint = testpoint - repmat(seedpoint, size(testpoint,1), 1);
[~, dists] = annsearch(pairwise', querypoint', K);
meandists = mean(dists, 1);
meandists = meandists';

end

function point = rect2point(rect, rule, data)
if rule == 1 % original
    point = (rect(:,1:3)+rect(:,4:6))/2;
elseif rule == 2 % normalized
    seedsize = data.normsize;
    point = (rect(:,1:3)+rect(:,4:6))/2./repmat(seedsize, size(rect,1), 1);
elseif rule == 3 % wall
    point = (rect(:,1:3)+rect(:,4:6))/2;
    pdim = rem(data.seedangle-1,2)+1;
    rdim = pdim + 3*(data.seedangle>=3);
    point(:,pdim) = rect(:,rdim);
end 
end