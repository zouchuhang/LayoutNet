function [ meandists ] = unaryHyperTest( testrects, unary, unarytype, data )
%UNARYHYPERTEST Summary of this function goes here
%   Detailed explanation goes here
if isempty(testrects)
    meandists = [];
    return;
end

testsize = testrects(:,4:6) - testrects(:,1:3);

if unarytype==2
    roomrect = data.roomrect;
    roomsize = roomrect(4:6) - roomrect(1:3);
    testsize = testsize./repmat(roomsize, size(testsize,1),1);
end

swap = data.testangle==2 | data.testangle==4;
t = testsize(swap,1);
testsize(swap,1) = testsize(swap,2);
testsize(swap,2) = t;

if size(unary,1)>5
    [~, dists] = annsearch(unary', testsize', 5);
    meandists = mean(dists, 1);
    meandists = meandists';
else
    dists = unary - repmat(testsize, size(unary,1), 1);
    dists = sqrt( sum(dists.^2,2));
    meandists = mean(dists);
end

% if unarytype==1
%     [~, dists] = annsearch(unary', testsize', 5);
%     meandists = mean(dists, 1);
%     meandists = meandists';
% elseif unarytype==2
%     
%     if data.testangle==2 || data.testangle==4
%         t = roomsize(1);
%         roomsize(1) = roomsize(2);
%         roomsize(2) = t;
%     end
%     testsize = testsize./repmat(roomsize, size(testsize,1),1);
%     [~, dists] = annsearch(unary', testsize', 5);
%     meandists = mean(dists, 1);
%     meandists = meandists';
% end

end

