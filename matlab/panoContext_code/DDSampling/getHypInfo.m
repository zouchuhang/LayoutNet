function [ comp ] = getHypInfo( hyps )
%GETHYPINFO Summary of this function goes here
%   Detailed explanation goes here
num_hyps = length(hyps);
comp = repmat(struct('objcenter',[],'objsize',[],'objtype',[]), num_hyps, 1 );
for i = 1:num_hyps
    if isempty(hyps{i})
        continue;
    end
    T = reshape(vertcat(hyps{i}.align)', 24, []);
    comp(i).objcenter = (T(1:3,:)+T(19:21,:))/2;
    comp(i).objcenter = comp(i).objcenter - ...
        repmat(comp(i).objcenter(:,1), 1, size(comp(i).objcenter, 2));
    comp(i).objsize = T(19:21,:)-T(1:3,:)+0.1;
    comp(i).objtype = horzcat(hyps{i}.objtype);
    comp(i).objaspectXY = comp(i).objsize(2,:)./comp(i).objsize(1,:);
    comp(i).objregionXY = comp(i).objsize(2,:).*comp(i).objsize(1,:);
    
    comp(i).objcenter = single(comp(i).objcenter);
    comp(i).objsize = single(comp(i).objsize);
    comp(i).objtype = single(comp(i).objtype);
    comp(i).objaspectXY = single(comp(i).objaspectXY);
    comp(i).objregionXY = single(comp(i).objregionXY);
end


end

