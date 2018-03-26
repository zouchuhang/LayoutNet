function roomhyp = sample_roomtype6(n, lines, vp, imgwidth, imgheight)
% sample n room hyps
% 3 ways to create hypothesis
% recipe1: 1,left 2,top, junc2type,1
% recipe2: 1,left 3,top,left, junc2type,2
% recipe3: 2,top 3,top,left, junc2type,3

lc = [lines.lineclass]; % 0,1,2,3
lr = [lines.leftorright]; % 1 if right, -1 if left, 0 if neither
tb = [lines.above_horizon]; % 1 if above, -1 if below, 0 if neither

oneleft = find((lc==1) & (lr==-1));
twotop = find((lc==2) & (tb==1));
threetopleft = find((lc==3) & (tb==1) & (lr==-1));

linecode{1} = oneleft;
linecode{2} = twotop;
linecode{3} = threetopleft;

% recipe: [linecode1 linecode2 junc2type]
recipe(1,:) = [1 2 1];
recipe(2,:) = [1 3 2];
recipe(3,:) = [2 3 3];

%
% for i = 1:length(linecode),
%     num_linecode(i) = length(linecode{i});
% end
if (length(linecode{1})<1 || length(linecode{2})<1) && ...
   (length(linecode{1})<1 || length(linecode{3})<1) && ...
   (length(linecode{2})<1 || length(linecode{3})<1)
    roomhyp = [];
    return; % fail
end

count = 0;
num_continue_without_progress = 0;
while count < n && num_continue_without_progress < 30
    use_recipe_num = randsample(size(recipe,1), 1); % 1, 2, or 3
    cur_recipe = recipe(use_recipe_num,:);

    l1 = linecode{cur_recipe(1)};
    l2 = linecode{cur_recipe(2)};
    if length(l1)<1 || length(l2)<1
        num_continue_without_progress = num_continue_without_progress + 1;
        continue;
    end

    ls1 = l1(randsample(length(l1),1));
    ls2 = l2(randsample(length(l2),1));
    [cornerpt degen] = line_intersect(...
                        lines(ls1).point1, lines(ls1).point2, ...
                        lines(ls2).point1, lines(ls2).point2);
    if degen==1
        num_continue_without_progress = num_continue_without_progress + 1;
        continue;
    end
    type = getjunctype(cornerpt, lines(ls1), lines(ls2));
    if type ~= cur_recipe(3)
        num_continue_without_progress = num_continue_without_progress + 1;
        continue;
    end

    % check if inside img
    MARGIN = 5;
    if ~is_in_image(cornerpt, imgwidth, imgheight, MARGIN)
        num_continue_without_progress = num_continue_without_progress + 1;
        continue;
    end
    
    % found a valid way to generate a room hyp
    count = count+1;
    num_continue_without_progress = 0;
    roomhyp(count).corner(4).pt = []; % up to 4 corners
    roomhyp(count).corner(1).pt = cornerpt;
    roomhyp(count).type = 6;
    
%     global img;
%     disp_vanish(img, lines([ls1 ls2]), vp);
%     plot(cornerpt(1), cornerpt(2), 'x');
%     pause;
%     close;
end

if ~exist('roomhyp','var')
    roomhyp = [];
end
