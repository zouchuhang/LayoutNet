function roomhyp = sample_roomtype5(n, lines, vp, imgwidth, imgheight)

% 6 ways to create hypotheses
% roomtype3 + 1,right
% roomtype3 + 3,top,right
% roomtype3 + 3,bot,right
% roomtype4 + 1,left
% roomtype4 + 3,top,left
% roomtype4 + 3,bot,left

lc = [lines.lineclass]; % 0,1,2,3
lr = [lines.leftorright]; % 1 if right, -1 if left, 0 if neither
tb = [lines.above_horizon]; % 1 if above, -1 if below, 0 if neither

oneright = find((lc==1) & (lr==1));
oneleft = find((lc==1) & (lr==-1));
threetopright = find((lc==3) & (tb==1) & (lr==1));
threebotright = find((lc==3) & (tb==-1) & (lr==1));
threetopleft = find((lc==3) & (tb==1) & (lr==-1));
threebotleft = find((lc==3) & (tb==-1) & (lr==-1));
linecode{1} = oneright;
linecode{2} = oneleft;
linecode{3} = threetopright;
linecode{4} = threebotright;
linecode{5} = threetopleft;
linecode{6} = threebotleft;

% recipe: [roomtypepart linecode oldcornerid newcornerid oldcornerid2 newcornerid2]
recipe(1,:) = [3 1 1 2 3 4]; % could be [3 1 3 4 1 2]
recipe(2,:) = [3 3 1 2 3 4];
recipe(3,:) = [3 4 3 4 1 2];
recipe(4,:) = [4 2 2 1 4 3]; % could be [4 2 4 3 2 1]
recipe(5,:) = [4 5 2 1 4 3];
recipe(6,:) = [4 6 4 3 2 1];

count = 0;
num_continue_without_progress = 0;
while count < n && num_continue_without_progress < 40
    % pick recipe
    cur_recipe = recipe(randsample(size(recipe,1), 1), :);

    % sample partial room
    if cur_recipe(1) == 3
        roomhyppart = sample_roomtype3(1, lines, vp, imgwidth, imgheight);
    elseif cur_recipe(1) == 4
        roomhyppart = sample_roomtype4(1, lines, vp, imgwidth, imgheight);
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
    
    % first new corner 
    oldpoint = roomhyppart.corner(cur_recipe(3)).pt;
    [newpoint1 degen] = line_intersect(vp{2},oldpoint,...
        lines(ls1).point1, lines(ls1).point2);
    if degen==1
        num_continue_without_progress = num_continue_without_progress + 1;
        continue;
    end
    % check if inside img
    MARGIN = 5;
    if ~is_in_image(newpoint1, imgwidth, imgheight, MARGIN)
        num_continue_without_progress = num_continue_without_progress + 1;
        continue;
    end
    
    % second new corner
    oldpoint2 = roomhyppart.corner(cur_recipe(5)).pt;
    [newpoint2 degen] = line_intersect(vp{2},oldpoint2, ...
        newpoint1, vp{1});
    if degen==1
        num_continue_without_progress = num_continue_without_progress + 1;
        continue;
    end
    % check if inside img
    MARGIN = 5;
    if ~is_in_image(newpoint2, imgwidth, imgheight, MARGIN)
        num_continue_without_progress = num_continue_without_progress + 1;
        continue;
    end

    % TODO: check location of new line sample
    %

    % add to hypothesis
    count = count + 1;
    num_continue_without_progress = 0;
    roomhyp(count) = roomhyppart;
    roomhyp(count).corner(cur_recipe(4)).pt = newpoint1;
    roomhyp(count).corner(cur_recipe(6)).pt = newpoint2;
    roomhyp(count).type = 5;
    
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
