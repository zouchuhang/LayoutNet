function vissurf = get_cube_visible_surface(cube, vp)

% cube surfaceid
% botid = 1; % 1: bottom
% frtid = 2; % 2: front
% lftid = 3; % 3: left
% rhtid = 4; % 4: right
% bckid = 5; % 5: back
% topid = 6; % 6: top

% cube junc3 id
%   5---6   1------2
%  /   /|   |      |
% 1---2 8   3------4    ....
% |   |/     \    /
% 3---4       7--8
%
% 1: front,top,left
% 2: front,top,right
% ...
% 8: back,bottom,right

%% determine visible surfaces

regionid1 = getregionid(cube.junc3(1).pt(1), cube.junc3(1).pt(2), vp);
regionid2 = getregionid(cube.junc3(2).pt(1), cube.junc3(2).pt(2), vp);
regionid3 = getregionid(cube.junc3(3).pt(1), cube.junc3(3).pt(2), vp);
regionid4 = getregionid(cube.junc3(4).pt(1), cube.junc3(4).pt(2), vp);

% vissurf = zeros(1,5); % [1234(front), 1256(top), 3478(bot), 1357(left), 2468(right)]
vissurf = zeros(1,6); % 3478(bot), 1234(front), 1357(left), 2468(right), 5678(back), 1256(top)

%bot
if (regionid3==1 || regionid3==2) && (regionid4==1 || regionid4==2)
    vissurf(1) = 1;
end

%front -- always visible
vissurf(2) = 1;

%left
if (regionid1==2 || regionid1==4) && (regionid3==2 || regionid3==4)
    vissurf(3) = 1;
end

%right
if (regionid2==1 || regionid2==3) || (regionid4==1 || regionid4==3)
    vissurf(4) = 1;
end

%back -- always invisible
vissurf(5) = 0;

%top
if (regionid1==3 || regionid1==4) && (regionid2==3 || regionid2==4)
    vissurf(6) = 1;
end


%% old
% % % regionid1 = getregionid(cube.junc3(1).pt(1), cube.junc3(1).pt(2), vp);
% % % regionid2 = getregionid(cube.junc3(2).pt(1), cube.junc3(2).pt(2), vp);
% % % regionid3 = getregionid(cube.junc3(3).pt(1), cube.junc3(3).pt(2), vp);
% % % regionid4 = getregionid(cube.junc3(4).pt(1), cube.junc3(4).pt(2), vp);
% % % 
% % % vissurf = zeros(1,5); % [1234(front), 1256(top), 3478(bot), 1357(left), 2468(right)]
% % % %front -- always visible
% % % vissurf(1) = 1;
% % % %top
% % % if (regionid1==3 || regionid1==4) && (regionid2==3 || regionid2==4)
% % %     vissurf(2) = 1;
% % % end
% % % %bot
% % % if (regionid3==1 || regionid3==2) && (regionid4==1 || regionid4==2)
% % %     vissurf(3) = 1;
% % % end
% % % %left
% % % if (regionid1==2 || regionid1==4) && (regionid3==2 || regionid3==4)
% % %     vissurf(4) = 1;
% % % end
% % % %right
% % % if (regionid2==1 || regionid2==3) || (regionid4==1 || regionid4==3)
% % %     vissurf(5) = 1;
% % % end
