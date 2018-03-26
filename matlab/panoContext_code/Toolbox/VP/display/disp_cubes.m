function h = disp_cubes(cubes, img, mode)

% cubes = drop_thin_cuboids(cubes, vp);

if nargin<3
    mode = 1;
end

if ~isempty(img)
    figure; imshow(img); hold on;
end

h = [];
for i = 1:length(cubes)
    h = [h disp_cube_2(cubes(i))];
    if mode==3
        pause; delete(h); h = [];
    elseif mode==2
        pause(0.01); delete(h); h = [];
    end
end

%%
function h = disp_cube_2(cube)

% color = rand(1,3);
color = [1 0.5 0];
% color = [1 1 0];
linewidth = 5;

h = zeros(1,12);
h(1) = plotline(cube.junc3(1).pt, cube.junc3(2).pt, color, linewidth);
h(2) = plotline(cube.junc3(1).pt, cube.junc3(3).pt, color, linewidth);
h(3) = plotline(cube.junc3(1).pt, cube.junc3(5).pt, color, linewidth);
h(4) = plotline(cube.junc3(2).pt, cube.junc3(4).pt, color, linewidth);
h(5) = plotline(cube.junc3(2).pt, cube.junc3(6).pt, color, linewidth);
h(6) = plotline(cube.junc3(3).pt, cube.junc3(4).pt, color, linewidth);
h(7) = plotline(cube.junc3(3).pt, cube.junc3(7).pt, color, linewidth);
h(8) = plotline(cube.junc3(4).pt, cube.junc3(8).pt, color, linewidth);
h(9) = plotline(cube.junc3(5).pt, cube.junc3(6).pt, color, linewidth);
h(10) = plotline(cube.junc3(5).pt, cube.junc3(7).pt, color, linewidth);
h(11) = plotline(cube.junc3(6).pt, cube.junc3(8).pt, color, linewidth);
h(12) = plotline(cube.junc3(7).pt, cube.junc3(8).pt, color, linewidth);

% for j = 1:8
%     h(12+j) = text(cube.junc3(j).pt(1), cube.junc3(j).pt(2), num2str(j), 'Color', [1 0 0]);
% end
