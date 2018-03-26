function roomhyp = sample_roomhyp(n, lines, vp, imgsize)

% room type:
%  1.      2.       3.        4.        5.           6.     7.
%                      \           /       \     /     \      /
%                       1--     --2         1---2       1-  -2
%    |        |         |         |         |   |       |    |
%    3--    --4         3--     --4         3---4             
%   /          \       /           \       /     \            
%

%%
roomhyp = [];
roomhyp = [roomhyp sample_roomtype1(n*2/10, lines, vp, imgsize(2), imgsize(1))];
roomhyp = [roomhyp sample_roomtype2(n*2/10, lines, vp, imgsize(2), imgsize(1))];

roomhyp = [roomhyp sample_roomtype6(n*2/10, lines, vp, imgsize(2), imgsize(1))];
roomhyp = [roomhyp sample_roomtype7(n*2/10, lines, vp, imgsize(2), imgsize(1))];

roomhyp = [roomhyp sample_roomtype3(n*2/10, lines, vp, imgsize(2), imgsize(1))];
roomhyp = [roomhyp sample_roomtype4(n*2/10, lines, vp, imgsize(2), imgsize(1))];

roomhyp = [roomhyp sample_roomtype5(n*4/10, lines, vp, imgsize(2), imgsize(1))];


%%
roomhyp = compute_box_from_room(roomhyp, vp, imgsize(2), imgsize(1));

%%
roomhyp = drop_roomhyp12_edge(roomhyp, imgsize(2));

