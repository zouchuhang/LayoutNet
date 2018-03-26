function [ selvalid ] = softSelect( hypScores, thres, soft )
%SOFTSELECT Summary of this function goes here
%   Detailed explanation goes here
nm = (thres-hypScores)./soft;
rn = rand(length(nm),1);
selvalid = rn>nm;

end

