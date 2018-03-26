function disp_lines(rectimg, lines, label, showimage)

if ~exist('showimage','var') || showimage==1
	figure; imshow(rectimg); hold on;
else
	hold on;
end
    
if nargin<3
    for k = 1:length(lines)
        xy = [lines(k).point1; lines(k).point2];
%         col = [0 0 1];
        col = rand(1,3);
        plot(xy(:,1),xy(:,2), 'LineWidth',2, 'Color', col);
        plot(xy(1,1), xy(1,2), 'x', 'LineWidth',2, 'Color', col);
        text(mean(xy(:,1)),mean(xy(:,2)), num2str(k), 'Color', col);
    end
else    
    for k = 1:length(lines)
        xy = [lines(k).point1; lines(k).point2];
text(xy(1,1), xy(1,2), sprintf('%d',k));
        if label(k) == 1
            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','red');
        elseif label(k) == -1
            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');
        elseif label(k) == 0
            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
        end
    end
end

