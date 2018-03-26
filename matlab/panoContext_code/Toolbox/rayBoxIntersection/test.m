% TEST:  Ray/box intersection using the Smits' algorithm

% origin    = [10.2, 10.2, 6.2];
% direction = [-0.5, -1, -1];
% vmin      = [-2.5, -4.5, -6.5];  % vertex min
% vmax      = [ 5,  5,  5];        % vertex max

%origin    = [0, 4, 2];
%direction = [0.213, -0.436, 0.873];
%vmin      = [-1 2 1];        % vertex min
%vmax      = [ 3 3 3];        % vertex max

origin = [0 0 0];
direction = [-1 -1 -1];
vmin = [2 2 2];
vmax = [4 4 4];

[flag ,tmin] = rayBoxIntersection(origin, direction, vmin, vmax);
intersection = origin + tmin*direction;

figure;
hold on;
grid on;

    % box (voxel)
    vertices = [vmax(1) vmin(2) vmin(3); vmax(1) vmax(2) vmin(3); vmin(1) vmax(2) vmin(3); vmin(1) vmax(2) vmax(3); vmin(1) vmin(2) vmax(3); vmax(1) vmin(2) vmax(3); vmin; vmax ];
    faces = [1 2 3 7; 1 2 8 6; 1 6 5 7; 7 5 4 3; 2 8 4 3; 8 6 5 4];
    h= patch('Vertices',vertices,'Faces',faces,'FaceColor','green');
%     set(h,'FaceAlpha',0.5);
  
    % origin
    text(origin(1), origin(2), origin(3), 'origin');
    plot3(origin(1), origin(2), origin(3), 'k.', 'MarkerSize', 10);

    % direction
    quiver3(origin(1), origin(2), origin(3), direction(1), direction(2), direction(3), 15);

    % intersection 
    plot3(intersection(1), intersection(2), intersection(3), 'r.', 'MarkerSize', 15);

    view(60,30);
    axis tight;
    xlabel('x');
    ylabel('y');
    zlabel('z');
