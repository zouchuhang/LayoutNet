function h = disp_room(roomhyp, img, showall, color, linewidth, linestyle)

% room type:
%  1.      2.       3.        4.        5.                    
%                      \           /       \     /            
%                       1--     --2         1---2             
%    |        |         |         |         |   |             
%    3--    --4         3--     --4         3---4             
%   /          \       /           \       /     \            
%

if ~exist('showall', 'var')
    showall = 1;
end
if ~exist('linestyle', 'var')
    linestyle = '-';
end
if ~exist('linewidth', 'var')
    linewidth = 5;
end

if ~isempty(img)
    figure; imshow(img); hold on;
end

for j = 1:length(roomhyp)
    if ~exist('color', 'var')
            color = [0 0 1];
    elseif showall==1
        if length(roomhyp)==1
            color = [0 0 1];
        else
            color = rand(1,3);
        end
    end

    
    h = [];
    if roomhyp(j).type==1
        pt = cat(1, roomhyp(j).box(1).p2, roomhyp(j).box(1).p4);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(2).p1, roomhyp(j).box(2).p2);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(2).p2, roomhyp(j).box(2).p4);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
    elseif roomhyp(j).type==2
        pt = cat(1, roomhyp(j).box(1).p2, roomhyp(j).box(1).p4);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(2).p1, roomhyp(j).box(2).p2);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(2).p2, roomhyp(j).box(2).p4);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
    elseif roomhyp(j).type==3
        pt = cat(1, roomhyp(j).box(1).p1, roomhyp(j).box(1).p3);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(1).p2, roomhyp(j).box(1).p4);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(2).p1, roomhyp(j).box(2).p2);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(2).p1, roomhyp(j).box(2).p3);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(2).p2, roomhyp(j).box(2).p4);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
    elseif roomhyp(j).type==4
        pt = cat(1, roomhyp(j).box(1).p1, roomhyp(j).box(1).p3);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(1).p2, roomhyp(j).box(1).p4);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(1).p3, roomhyp(j).box(1).p4);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(2).p1, roomhyp(j).box(2).p3);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(2).p2, roomhyp(j).box(2).p4);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
    elseif roomhyp(j).type==5
        pt = cat(1, roomhyp(j).box(1).p1, roomhyp(j).box(1).p3);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(1).p2, roomhyp(j).box(1).p4);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(2).p1, roomhyp(j).box(2).p2);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(2).p1, roomhyp(j).box(2).p3);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(2).p2, roomhyp(j).box(2).p4);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(2).p3, roomhyp(j).box(2).p4);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(3).p1, roomhyp(j).box(3).p3);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(3).p2, roomhyp(j).box(3).p4);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
    elseif roomhyp(j).type==6
        pt = cat(1, roomhyp(j).box(1).p1, roomhyp(j).box(1).p3);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(1).p3, roomhyp(j).box(1).p4);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(2).p1, roomhyp(j).box(2).p3);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
    elseif roomhyp(j).type==7
        pt = cat(1, roomhyp(j).box(1).p1, roomhyp(j).box(1).p3);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(1).p3, roomhyp(j).box(1).p4);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
        pt = cat(1, roomhyp(j).box(2).p1, roomhyp(j).box(2).p3);
        h = [h plot(pt(:,1), pt(:,2), linestyle, 'LineWidth', linewidth, 'Color', color)];
    end
    
% % % %     if roomhyp(j).type==3 || roomhyp(j).type==5
% % % %         pt = cat(1, roomhyp(j).box(1).p1, roomhyp(j).box(1).p3);
% % % %         h = [h plot(pt(:,1), pt(:,2), 'LineWidth', 2, 'Color', color)];
% % % %     end
% % % %     if roomhyp(j).type==1 || roomhyp(j).type==3 || roomhyp(j).type==5
% % % %         pt = cat(1, roomhyp(j).box(1).p2, roomhyp(j).box(1).p4);
% % % %         h = [h plot(pt(:,1), pt(:,2), 'LineWidth', 2, 'Color', color)];
% % % %     end
% % % %     if roomhyp(j).type==1 || roomhyp(j).type==3 || roomhyp(j).type==5
% % % %         pt = cat(1, roomhyp(j).box(2).p1, roomhyp(j).box(2).p2);
% % % %         h = [h plot(pt(:,1), pt(:,2), 'LineWidth', 2, 'Color', color)];
% % % %     end
% % % %     pt = cat(1, roomhyp(j).box(2).p1, roomhyp(j).box(2).p3);
% % % %     h = [h plot(pt(:,1), pt(:,2), 'LineWidth', 2, 'Color', color)];
% % % %     pt = cat(1, roomhyp(j).box(2).p2, roomhyp(j).box(2).p4);
% % % %     h = [h plot(pt(:,1), pt(:,2), 'LineWidth', 2, 'Color', color)];
% % % %     if roomhyp(j).type==2 || roomhyp(j).type==4 || roomhyp(j).type==5
% % % %         pt = cat(1, roomhyp(j).box(2).p3, roomhyp(j).box(2).p4);
% % % %         h = [h plot(pt(:,1), pt(:,2), 'LineWidth', 2, 'Color', color)];
% % % %     end
% % % %     if roomhyp(j).type==4 || roomhyp(j).type==5
% % % %         pt = cat(1, roomhyp(j).box(3).p1, roomhyp(j).box(3).p3);
% % % %         h = [h plot(pt(:,1), pt(:,2), 'LineWidth', 2, 'Color', color)];
% % % %     end
% % % %     if roomhyp(j).type==2 || roomhyp(j).type==4 || roomhyp(j).type==5
% % % %         pt = cat(1, roomhyp(j).box(3).p2, roomhyp(j).box(3).p4);
% % % %         h = [h plot(pt(:,1), pt(:,2), 'LineWidth', 2, 'Color', color)];
% % % %     end

    if showall == 2
        pause(0.01); delete(h);
    elseif showall == 3
        pause; delete(h);
    end
end


% % % % for j = 1:length(roomhyp)
% % % %     h = [];
% % % % %     for i = 1:4
% % % % %         if ~isempty(roomhyp(j).corner(i).pt)
% % % % %             hh = text(roomhyp(j).corner(i).pt(1), roomhyp(j).corner(i).pt(2), num2str(i), 'FontSize', 20);
% % % % %             h = [h hh];
% % % % %             pt = cat(1,roomhyp(j).corner([1 3]).pt);
% % % %         if roomhyp(j).type==1
% % % %             cid = 3;
% % % %             h = text(roomhyp(j).corner(cid).pt(1), roomhyp(j).corner(cid).pt(2),...
% % % %                 num2str(roomhyp(j).type), 'FontSize', 20);
% % % %         elseif roomhyp(j).type==2
% % % %             cid = 4;
% % % %             h = text(roomhyp(j).corner(cid).pt(1), roomhyp(j).corner(cid).pt(2),...
% % % %                 num2str(roomhyp(j).type), 'FontSize', 20);
% % % %         elseif roomhyp(j).type==3
% % % %             pt = cat(1,roomhyp(j).corner([1 3]).pt);
% % % %             h = plot(pt(:,1), pt(:,2), '-x','LineWidth', 2, 'MarkerSize', 5);
% % % %         elseif roomhyp(j).type==4
% % % %             pt = cat(1,roomhyp(j).corner([2 4]).pt);
% % % %             h = plot(pt(:,1), pt(:,2), '-x','LineWidth', 2, 'MarkerSize', 5);
% % % %         elseif roomhyp(j).type==5
% % % %             pt = cat(1,roomhyp(j).corner([1 2 4 3 1]).pt);
% % % %             h = plot(pt(:,1), pt(:,2), '-x','LineWidth', 2, 'MarkerSize', 5);
% % % %         end
% % % % 
% % % % %         end
% % % % %     end
% % % % %     pause;
% % % % %     delete(h);
% % % % end
