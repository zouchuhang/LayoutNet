function [ OUT_HYPS, DEL_TYPE, DEL_OBID ] = deleteObject( objscore, ALL_HYPS )
%DELETEOBJECT Summary of this function goes here
%   Detailed explanation goes here
numhyps = length(ALL_HYPS);
OUT_HYPS = ALL_HYPS;
DEL_TYPE = zeros(numhyps,1);
DEL_OBID = zeros(numhyps,1);
for i = 1:numhyps
%     fprintf('%d\n',i);
    sel_hyps = ALL_HYPS(i).sel_hyps;
    objIDs = zeros(0,1);
    typIDs = zeros(0,1);
    objScr = zeros(0,1);
    for j = 1:length(sel_hyps)
        if sel_hyps(j).fixed
            objIDs = [objIDs; sel_hyps(j).selID'];
            typIDs = [typIDs; j*ones(length(sel_hyps(j).selID),1)];
            objScr = [objScr; objscore(sel_hyps(j).selID, j)];
        end
    end
    [~,I] = sort(objScr,'descend');
    delID = I(randp([1:length(I)], 1, 1));
    sel_hyps(typIDs(delID)).selID = setdiff(sel_hyps(typIDs(delID)).selID, objIDs(delID));
    if isempty(sel_hyps(typIDs(delID)).selID)
        sel_hyps(typIDs(delID)).fixed = false;
    end
    DEL_TYPE(i) = typIDs(delID);
    DEL_OBID(i) = objIDs(delID);
    OUT_HYPS(i).sel_hyps = sel_hyps;
end

end

