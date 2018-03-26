function [ MINSCORE, MINTRANS] = compRoomMatchScore( HYPS, VALID_TRANS_GNDS )
%COMPROOMMATCHSCORE Summary of this function goes here
%   Detailed explanation goes here
HYP_COMP = getHypInfo(HYPS);
MINSCORE = zeros(length(HYPS), length(VALID_TRANS_GNDS));
MINTRANS = zeros(length(HYPS), length(VALID_TRANS_GNDS));
for j = 1:length(VALID_TRANS_GNDS)
%     fprintf('GND: %d/%d\n', j, length(VALID_TRANS_GNDS));
    fprintf('>>%d', j);
    tCOST = roomAlignmentMex(HYP_COMP, VALID_TRANS_GNDS(j).SAMPLE_COMP, 3, 0);
    [B,I] = min( tCOST, [], 2);
    MINSCORE(:,j) = B;
    MINTRANS(:,j) = I;
end

end

