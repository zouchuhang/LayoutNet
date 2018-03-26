function roomhyp_drop = drop_roomhyp12_edge(roomhyp, imgwidth)
% Deletes roomhyp type 1 and 2 where the top of the wall boundary touches
% a side of the image.
% This is a temporary fix. get_labelimg_frombox.m did not produce correct
% results with such case...
%
% e.g. +----------------+
%      |                |
%      |               /|  <-- here.. the vertical edge ends on image side
%      |              / |
%      |--------------\ |
%      |               \|
%      +----------------+

MARGIN = 3;

idxtodrop = zeros(1, length(roomhyp));
for i = 1:length(roomhyp)
    if roomhyp(i).type==1
        if roomhyp(i).box(1).p3(1) < 1+MARGIN
            idxtodrop(i) = 1;
        end
    elseif roomhyp(i).type==2
        if roomhyp(i).box(2).p1(1) > imgwidth-MARGIN
            idxtodrop(i) = 1;
        end
    end
end

roomhyp_drop = roomhyp(~idxtodrop);

