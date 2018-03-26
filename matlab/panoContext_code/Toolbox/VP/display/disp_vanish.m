function h = disp_vanish(rectimg, lines, vp)
% disp_vanish.m

%%
h = figure; imshow(rectimg); hold on;

%%
if length(vp)==3 && ~isempty(vp{1}) && ~isempty(vp{2}) && ~isempty(vp{3})
    plot([vp{1}(1) vp{2}(1) vp{3}(1)], [vp{1}(2) vp{2}(2) vp{3}(2)], 'bx', 'Markersize', 20, 'LineWidth', 4);
    % plot([vp{2}(1) vp{3}(1)], [vp{2}(2) vp{3}(2)], 'bx', 'Markersize', 20, 'LineWidth', 4);
    try
        [p1ext p2ext] = extline(vp{1}, vp{3}, size(rectimg,2), size(rectimg,1));
        % plot([p1ext(1) p2ext(1)], [p1ext(2) p2ext(2)], 'w', 'LineWidth', 3);
        plot([p1ext(1) p2ext(1)], [p1ext(2) p2ext(2)], 'k', 'LineWidth', 3);
    catch
        % extline might fail when both vp's are outside of image
    end
    try
        [p1ext p2ext] = extline(vp{2}, vp{3}, size(rectimg,2), size(rectimg,1));
        % plot([p1ext(1) p2ext(1)], [p1ext(2) p2ext(2)], 'w', 'LineWidth', 3);
        plot([p1ext(1) p2ext(1)], [p1ext(2) p2ext(2)], 'k', 'LineWidth', 3);
    catch
        % extline might fail when both vp's are outside of image
    end
end

%% undistorted
% figure; hold on;
% plot(vp1(1), vp1(2), 'r.', 'MarkerSize',20);
% plot(vp2(1), vp2(2), 'g.', 'MarkerSize',20);
% plot(vp3(1), vp3(2), 'b.', 'MarkerSize',20);
for k = 1:length(lines)
% 	xy = [lines(k).point1_dist; lines(k).point2_dist];
	xy = [lines(k).point1; lines(k).point2];
	if lines(k).lineclass == 1
		plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','red');
	elseif lines(k).lineclass == 2
		plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
	elseif lines(k).lineclass == 3
		plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');
	else
		plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','yellow');
	end
end
% disp(round([vp1; vp2; vp3]));
% disp([normal_vec{1}; normal_vec{2}; normal_vec{3}]);

%% distorted
% figure; imshow(rectimg); hold on;
% plot(vp1(1), vp1(2), 'r.', 'MarkerSize',20);
% plot(vp2(1), vp2(2), 'g.', 'MarkerSize',20);
% plot(vp3(1), vp3(2), 'b.', 'MarkerSize',20);
% for k = 1:length(lines)
% 	xy = [lines(k).point1; lines(k).point2];
% 	if lines(k).lineclass1 == 1
% 		plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','red');
% 	elseif lines(k).lineclass2 == 1
% 		plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% 	elseif lines(k).lineclass3 == 1
% 		plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');
% 	else
% 		plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','black');
% 	end
% end
% disp(round([vp1; vp2; vp3]));
% disp([n1; n2; n3]);
% 
