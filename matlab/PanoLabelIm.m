function PanoLabelIm(id, ANNO_ALL, config, options)
% load image
im_name = config.IMGLIST(id).name;
im = imread([config.folderName im_name '/' im_name '.jpg']);
[im_h, im_w, im_dim] = size(im);
im_o = im2double(im);
% alignment
if 1
imgSize = 320;
qError = 0.7;
global config;
config.NewProjFolder = './panoContext_code/';
config.lsdLocation = [config.NewProjFolder 'VpEstimation/lsd '];
config.bufferFolder = './Buffer_gt/';
if ~exist(config.bufferFolder, 'dir')
    mkdir(config.bufferFolder);
end
[ olines, vp, views, edges, panoEdge, score, angle] = panoEdgeDetection( im_o, imgSize, qError);
save(['./label_vp/pc/' im_name '.mat'], 'vp');
%keyboard
return
end
load(['./label_vp/pc/' im_name '.mat']);
[rotImg, R] = rotatePanorama(im_o, vp(3:-1:1,:));

if exist(['./label_cor/pc/' im_name '.mat'], 'file')
    gt = load(['./label_cor/pc/' im_name '.mat']);
    gt = gt.cor;
    keyboard
else

r_id = find([ANNO_ALL(id).ANNO3D.objects.type] == 3);
    pt = ANNO_ALL(id).ANNO3D.objects(r_id).points *ANNO_ALL(id).ANNO3D.R;
    pt = pt*R'; % algin pt also
    [ uv ] = xyz2uvN( pt ); % transform xyz to im coordinates
    [ uvcoord ] = uv2coords( uv, im_w, im_h );
%sort
[~,gt_id] = sort(uvcoord(:,1));
gt = uvcoord(gt_id,:);
for i = 1:2:size(gt,1)
    gt_t = gt(i:i+1,:);
    [~, gt_t_id] = sort(gt_t(:,2));
    gt(i:i+1,:) = gt_t(gt_t_id,:);
end
end
gt_all = [gt;
               gt(1,:);gt(3,:);gt(3,:);gt(5,:);
               gt(5,:);gt(7,:);gt(7,:);gt(1,:);
               gt(2,:);gt(4,:);gt(4,:);gt(6,:);
               gt(6,:);gt(8,:);gt(8,:);gt(2,:)];
    [ uv ] = coords2uv( gt_all, im_w, im_h );
    [ xyz ] = uv2xyzN( uv);
    [ lines ] = lineFromTwoPoint( xyz(1:2:end,:), xyz(2:2:end,:) );
    [ panoEdgeC2 ] = paintParameterLine2(lines, im_w, im_h, rotImg, 2 );

figure(1), axis off, hold off,
ha = tight_subplot(1, 1);
axes(ha(1)), imagesc(panoEdgeC2), axis image, hold on

if 1 % corner debug
disp('Click y if label this image');
[~,~,b] = ginput(1);
if b ~= 'y'
    return
end

% Allow user to input line segments; compute centers, directions, lengths
%disp('Draw Pano lines...')

% find concave lines
lines = [];
while 1
    disp(' ')
    disp('Draw concave lines or q to stop')
    [x1,y1,b] = ginput(1);    
    if b=='q'        
        break;
    end
    [x2,y2, ~] = ginput(1);
    % single line fitting
    vert = (x1+x2)/2;
    axes(ha(1)); h1 = plot([vert vert], [1 im_h], 'b');
    disp('Click y if accept this line')
    [~,~, b] = ginput(1);
    % save
    if b == 'y'
        lines(end+1,:) = [vert y1 1]; % 1 for concave signal
        lines(end+1,:) = [vert y2 1]; 
    else
        axes(ha(1)); delete(h1);
    end
end

% find convex lines
while 1
    disp(' ')
    disp('Draw convex lines or q to stop')
    [x1,y1,b] = ginput(1);    
    if b=='q'        
        break;
    end
    [x2,y2, ~] = ginput(1);
    % single line fitting
    vert = (x1+x2)/2;
    axes(ha(1)); h1 = plot([vert vert], [1 im_h], 'r');
    disp('Click y if accept this line')
    [~,~, b] = ginput(1);
    % save
    if b == 'y'
        lines(end+1,:) = [vert y1 2]; % 2 for convex signal
        lines(end+1,:) = [vert y2 2]; 
    else
        axes(ha(1)); delete(h1);
    end
end


% sort lines
[~,l_id] = sort(lines(:,1));
lines = lines(l_id,:);
for i = 2:2:size(lines,1)
    sub_l = lines(i-1:i,:);
    [~,l_id] = sort(sub_l,1);
    lines(i-1,:) = sub_l(l_id(1),:);
    lines(i,:) = sub_l(l_id(2),:);
end
end

% solve for pano
if size(lines,1) == 8
    [wall_d, x, f] = pano_line_solver(lines, im_w, options); % bos-shape case
else if size(lines,1) == 12
       [wall_d, x, f] = pano_line_solver_6(lines, im_w, options); % L-shape case 
    else if size(lines,1) == 16
            [wall_d, x, f] = pano_line_solver_8(lines, im_w, options); % T-shape case 
        end
    end
end

lines_ = lines;
% annotate box height
cor = [];
%c_h = 1.7;
flag = 0;
while 1
    disp('Draw corners or q to stop')
    [x1,y1,b] = ginput(1); 
    axes(ha(1)); h1 = plot(x1, y1, 'xr');   
    if b=='q'        
        break;
    end
    [x2,y2] = ginput(1);
    axes(ha(1)); h3 = plot(x2, y2, 'xr');
    disp('Click y if accept this corner');
    [~,~,b] = ginput(1);
    if b == 'y'
        cor(end+1,:) = [x1 y1];
        cor(end+1,:) = [x2 y2]; 
        % find closet wall
        %keyboard
        vert = (x1+x2)/2;
        cor_d = abs(lines(:,1) - vert);
        [~, l_id] = min(cor_d);
        cor = [lines(l_id,1) cor(1,2); lines(l_id,1) cor(2,2)];
        % solve for all walls
        [~, cor_id] = sort(cor(:,2));
        cor = cor(cor_id,:);
        theta_y = pi*cor(2,2)/im_h-pi/2;
        d = abs(cot(theta_y));lines_ = lines;%c_h*abs(cot(theta_y));
        theta_y_ = pi/2 - pi*cor(1,2)/im_h;
        h = d*tan(theta_y_);%c_h + d*tan(theta_y_);
        d_m = d/wall_d((l_id+1)/2);
        line_ya = [];
        line_ya_ = [];
        for i = 1:2:size(lines,1)
            if i == l_id
                continue
            end
            line_d = d_m*wall_d((i+1)/2);
            line_theta_y = acot(line_d);%acot(line_d/c_h);
            line_y = (line_theta_y + pi/2)*im_h/pi;
            cor = [cor;lines(i,1) line_y];
            line_y_ = atan(h/line_d);%atan((h-c_h)/line_d);
            line_y_ = (pi/2-line_y_)*im_h/pi;
            cor = [cor;lines(i,1) line_y_ ];
            axes(ha(1)); h1 = plot([lines(i,1) lines(i,1)], [line_y line_y_], 'xr');
        end
        % sort cor
        [~, cor_id] = sort(cor(:,1));
        cor = cor(cor_id,:);
        for i = 1:2:size(cor,1)
                cor_sub = cor(i:i+1,:);
                [~, cor_id] = sort(cor_sub(:,2));
                cor(i,:) = cor_sub(cor_id(1),:);
                cor(i+1,:) = cor_sub(cor_id(2),:);
        end

        % visual full
        if size(cor,1) == 8
            cor_all = [cor;
                cor(1,:);cor(3,:);cor(3,:);cor(5,:);
                cor(5,:);cor(7,:);cor(7,:);cor(1,:);
                cor(2,:);cor(4,:);cor(4,:);cor(6,:);
                cor(6,:);cor(8,:);cor(8,:);cor(2,:)];
        else if size(cor,1) == 12
            cor_all = [cor;
                cor(1,:);cor(3,:);cor(3,:);cor(5,:);
                cor(5,:);cor(7,:);cor(7,:);cor(9,:);
                cor(9,:);cor(11,:);cor(11,:);cor(1,:);
                cor(2,:);cor(4,:);cor(4,:);cor(6,:);
                cor(6,:);cor(8,:);cor(8,:);cor(10,:);
                cor(10,:);cor(12,:);cor(12,:);cor(2,:)];
		else if size(cor,1) == 16
		    cor_all = [cor;
                	cor(1,:);cor(3,:);cor(3,:);cor(5,:);
                	cor(5,:);cor(7,:);cor(7,:);cor(9,:);
                	cor(9,:);cor(11,:);cor(11,:);cor(13,:);
                        cor(13,:);cor(15,:);cor(15,:);cor(1,:);
                	cor(2,:);cor(4,:);cor(4,:);cor(6,:);
                	cor(6,:);cor(8,:);cor(8,:);cor(10,:);
                	cor(10,:);cor(12,:);cor(12,:);cor(14,:);
			cor(14,:);cor(16,:);cor(16,:);cor(2,:)];
		end
            end
        end
        [ uv ] = coords2uv( cor_all, im_w, im_h );
        [ xyz ] = uv2xyzN( uv);
        [ lines ] = lineFromTwoPoint( xyz(1:2:end,:), xyz(2:2:end,:) );
        [ panoEdgeC ] = paintParameterLine2(lines, im_w, im_h, rotImg, 1 );
        axes(ha(1));imagesc(panoEdgeC); axis image; hold on;
        for i = 1:size(line,1)
            axes(ha(1)); plot([lines(i,1) lines(i,1)], [1 im_h], 'b');
            axes(ha(1)); plot([lines(i,1) lines(i,1)], [line_y line_y_], 'xr');
        end
        
        disp('Click y if accept this box');
        [~,~,b] = ginput(1);
        if b == 'y'
            flag = 1;
            break
        else
            axes(ha(1)); imagesc(rotImg), axis image, hold on 
            for i = 1:size(lines,1)
                axes(ha(1)); plot([lines(i,1) lines(i,1)], [1 im_h], 'b');
            end
            cor = [];
        end
    else
        axes(ha(1)); delete(h1); delete(h3);
    end
    if flag
        break
    end
end

if 0
% sort cor
[~, cor_id] = sort(cor(:,1));
cor = cor(cor_id,:);
for i = 1:2:size(cor,1)
    cor_sub = cor(i:i+1,:);
    [~, cor_id] = sort(cor_sub(:,2));
    cor(i,:) = cor_sub(cor_id(1),:);
    cor(i+1,:) = cor_sub(cor_id(2),:);
end

% visual full
cor_all = [cor;
               cor(1,:);cor(3,:);cor(3,:);cor(5,:);
               cor(5,:);cor(7,:);cor(7,:);cor(1,:);
               cor(2,:);cor(4,:);cor(4,:);cor(6,:);
               cor(6,:);cor(8,:);cor(8,:);cor(2,:)];
[ uv ] = coords2uv( cor_all, im_w, im_h );
[ xyz ] = uv2xyzN( uv);
[ lines ] = lineFromTwoPoint( xyz(1:2:end,:), xyz(2:2:end,:) );
[ panoEdgeC ] = paintParameterLine2(lines, im_w, im_h, rotImg, 1 );
axes(ha(1));imagesc(panoEdgeC); axis image;
end

cor = [cor lines_(:,3)];

imwrite(panoEdgeC, ['./label_vis/pc/' im_name '.png'], 'png');
save(['./label_cor/pc/' im_name '.mat'], 'cor');
