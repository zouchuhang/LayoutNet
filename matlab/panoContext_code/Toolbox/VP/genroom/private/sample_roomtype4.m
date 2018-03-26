function roomhyp = sample_roomtype4(n, lines, vp, imgwidth, imgheight)

% 4 ways to create hypotheses
% roomtype2 + 3,top,right
% roomtype2 + 2,top
% roomtype7 + 3,bot,right
% roomtype7 + 2,bot

lc = [lines.lineclass]; % 0,1,2,3
lr = [lines.leftorright]; % 1 if right, -1 if left, 0 if neither
tb = [lines.above_horizon]; % 1 if above, -1 if below, 0 if neither

twotop = find((lc==2) & (tb==1));
twobot = find((lc==2) & (tb==-1));
threetopright = find((lc==3) & (tb==1) & (lr==1));
threebotright = find((lc==3) & (tb==-1) & (lr==1));

linecode{1} = twotop;
linecode{2} = twobot;
linecode{3} = threetopright;
linecode{4} = threebotright;

% recipe: [roomtypepart linecode oldcornerid newcornerid]
recipe(1,:) = [2 3 4 2];
recipe(2,:) = [2 1 4 2];
recipe(3,:) = [7 4 2 4];
recipe(4,:) = [7 2 2 4];

% roomhyp = [];
count = 0;
num_continue_without_progress = 0;
while count < n && num_continue_without_progress < 30
% pick recipe
    cur_recipe = recipe(randsample(size(recipe,1), 1), :);
    
    % sample partial room
    if cur_recipe(1) == 2
        roomhyppart = sample_roomtype2(1, lines, vp, imgwidth, imgheight);
    elseif cur_recipe(1) == 7
        roomhyppart = sample_roomtype7(1, lines, vp, imgwidth, imgheight);
    end
    
    if isempty(roomhyppart)
        num_continue_without_progress = num_continue_without_progress + 1;
        continue;
    end

    % sample lines
    l1 = linecode{cur_recipe(2)};
    if length(l1)<1
        num_continue_without_progress = num_continue_without_progress + 1;
        continue;
    end
    ls1 = l1(randsample(length(l1), 1));
    
    % new corner
    oldpoint = roomhyppart.corner(cur_recipe(3)).pt;
    [newpoint degen] = line_intersect(vp{1},oldpoint,...
        lines(ls1).point1, lines(ls1).point2);
    if degen==1
        num_continue_without_progress = num_continue_without_progress + 1;
        continue;
    end
    % check if inside img
    MARGIN = 5;
    if ~is_in_image(newpoint, imgwidth, imgheight, MARGIN)
        num_continue_without_progress = num_continue_without_progress + 1;
        continue;
    end
    
    % TODO: check location of new line sample
    %
    
    % add to hypothesis
    count = count + 1;
    num_continue_without_progress = 0;
    roomhyp(count) = roomhyppart;
    roomhyp(count).corner(cur_recipe(4)).pt = newpoint;
    roomhyp(count).type = 4;
    
%     global img;
%     disp_vanish(img, lines(ls1), vp);
%     plot(oldpoint(1), oldpoint(2), 'rx', 'MarkerSize',10, 'LineWidth',2);
%     plot(newpoint(1), newpoint(2), 'bx', 'MarkerSize',10, 'LineWidth',2);
%     pause;
%     close;
end

if ~exist('roomhyp','var')
    roomhyp = [];
end
