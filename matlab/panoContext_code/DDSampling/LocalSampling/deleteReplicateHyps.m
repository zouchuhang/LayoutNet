function [ valid ] = deleteReplicateHyps( hyps )
%DELETEREPLICATEHYPS Summary of this function goes here
%   Detailed explanation goes here
valid = true(length(hyps),1);
for i = 1:length(hyps)
    if ~valid(i)
        continue;
    end
    for j = i+1:length(hyps)
        if ~valid(j)
            continue;
        end
        sel_hyps1 = hyps(i).sel_hyps;
        sel_hyps2 = hyps(j).sel_hyps;
        v = true;
        for k = 1:12
            id1 = sort(sel_hyps1(k).selID);
            id2 = sort(sel_hyps2(k).selID);
            if length(id1) ~= length(id2)
                v = false; break;
            elseif ~isempty(id1) && ~isempty(id2)
                if sum(abs(id1-id2))~=0
                    v = false; break;
                end
            end
            
        end
        if v
            valid(j) = false;
        end
    end
end

end

