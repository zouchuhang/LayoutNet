% sample script to preprocess the panorama, edge map, corner map and box parameters

global config;
% train & test split
% NOTE: download PanoContext data and change the below path
config.folderName = '../panoContext_data/bedroom/';
config.annotationfile = [config.folderName 'ANNO_ALL.mat'];
load(config.annotationfile);
load([config.folderName 'IMGLIST.mat']); % image name
config.IMGLIST = IMGLIST;   
clear IMGLIST
anno_valid = false( 1, length(ANNO_ALL));
for i = 1:length(anno_valid)
    anno_valid(i) = ~isempty(ANNO_ALL(i).ANNO3D) && ANNO_ALL(i).ANNO3D.b_singleroom;
end
%clear ANNO_ALL;
config.anno_valid = anno_valid;

config.valid_anno_id = find(anno_valid);
config.valid_data_id = 1:sum(anno_valid);

config.NewProjFolder = './panoContext_code/';
config.lsdLocation = [config.NewProjFolder 'VpEstimation/lsd '];
config.bufferFolder = './Buffer/';
if ~exist(config.bufferFolder, 'dir')
    mkdir(config.bufferFolder);
end

clear anno_valid
clear i

config.TEST_DATA_IDS = [4 5 8 14 16 18 22 25 28 33 36 38 44 48 68 72 99 111 112 118 ...
    122 124 137 142 159 196 199 201 205 211 215 222 224 231 238 255 ...
    256 264 274 275 286 315 328 333 341 342 354 363 368 370 401 404 417];
config.TRAIN_DATA_IDS = setdiff(config.valid_anno_id, config.TEST_DATA_IDS);

% params
h = fspecial('gaussian',141 ,20);
im_h = 512;
im_w = 1024;
visual = 0;


se1 = strel('line',6,0);
se2 = strel('line',6,90);

% get gt
for i = 1:numel(config.TRAIN_DATA_IDS)
    tic
    id = config.TRAIN_DATA_IDS(i);
    
    % read pano
    im_name = config.IMGLIST(id).name;
    disp(im_name)
    panoImg_ori = imread([config.folderName im_name '/' im_name '.jpg']);
    panoImg_ori = im2double(panoImg_ori);
    [height, width, ~] = size(panoImg_ori);
        
    load(['../data/pano_edge_tr_1024/vp/' im_name '.mat']);
    
    Img_small = imresize(panoImg_ori, [1024 2048]);
    [rotImg, R] = rotatePanorama(Img_small, vp(3:-1:1,:));
    panoImg = imresize(rotImg, [im_h im_w]);

    % no alignment
    %panoImg = imresize(panoImg_ori, [im_h im_w]);
    
    % read gt
    coords = load(['./label_cor/pc/' im_name '.mat']);
    coords = coords.cor;
    [ uv ] = coords2uv( coords, width, height );
    [ uvcoord ] = uv2coords( uv, im_w, im_h );

    [ pt ] = uv2xyzN( uv);
    pt = pt*R; % comment this line if no alignment
    [ uv ] = xyz2uvN( pt );
    [ uvcoord ] = uv2coords( uv, im_w, im_h );
    
    %keyboard
    kpmap = zeros(im_h,im_w,3);
    
    kpcol = zeros(1,im_w);
    % generate gt map
    [~,uvid] = sort(uvcoord(:,1));
    uvcoord = uvcoord(uvid,:);
    
    cor_id = uvcoord(1:8,:);
    for j = 1:2:size(cor_id(:,1))
        cor_tmp = cor_id(j:j+1,:);
        [~,cor_t_id] = sort(cor_tmp(:,2));
        cor_id(j,:) = cor_tmp(cor_t_id(1),:);
        cor_id(j+1,:) = cor_tmp(cor_t_id(2),:);
    end
    % record rotation res
    box_n_all = zeros(im_w, 12);
    cor_id_o = cor_id;
    for rot = 0:im_w-1
        cor_id = cor_id_o;
        cor_id(:,1) = cor_id_o(:,1) + rot;
        msk = cor_id(:,1) > im_w;
        cor_id(msk,1) = cor_id(msk,1) - im_w;
        if sum(msk>0)
           [~,corid] = sort(cor_id(:,1));
           cor_id = cor_id(corid,:);
        end
        % gen box gt
        cor_a_x = cor_id(2,1) - im_w/2;
        cor_a_y = (im_h - cor_id(2,2))-im_h/2;
        cor_a_y_ = (im_h - cor_id(1,2))-im_h/2;
        theta_x = 2*pi*cor_a_x/im_w;
        theta_y = pi*cor_a_y/im_h;
        theta_y_ = pi*cor_a_y_/im_h;
        r = abs(cot(theta_y)); %*ch
        cor_a_X = r*cos(-theta_x+pi/2);% *ch
        cor_a_Z = r*sin(-theta_x+pi/2);% *ch 
        cor_a_Y_ = r*tan(theta_y_);% +ch

        cor_b_x = cor_id(4,1) - im_w/2;
        cor_b_y = (im_h - cor_id(4,2))-im_h/2;
        cor_b_y_ = (im_h - cor_id(3,2))-im_h/2;
        theta_x = 2*pi*cor_b_x/im_w;
        theta_y = pi*cor_b_y/im_h;
        theta_y_ = pi*cor_b_y_/im_h;
        r = abs(cot(theta_y)); %*ch
        cor_b_X = r*cos(-theta_x+pi/2);% *ch
        cor_b_Z = r*sin(-theta_x+pi/2);% *ch 
        cor_b_Y_ = r*tan(theta_y_);% +ch

        cor_c_x = cor_id(6,1) - im_w/2;
        cor_c_y = (im_h - cor_id(6,2))-im_h/2;
        cor_c_y_ = (im_h - cor_id(5,2))-im_h/2;
        theta_x = 2*pi*cor_c_x/im_w;
        theta_y = pi*cor_c_y/im_h;
        theta_y_ = pi*cor_c_y_/im_h;
        r = abs(cot(theta_y)); %*ch
        cor_c_X = r*cos(-theta_x+pi/2);% *ch
        cor_c_Z = r*sin(-theta_x+pi/2);% *ch 
        cor_c_Y_ = r*tan(theta_y_);% +ch

        cor_d_x = cor_id(8,1) - im_w/2;
        cor_d_y = (im_h - cor_id(8,2))-im_h/2;
        cor_d_y_ = (im_h - cor_id(7,2))-im_h/2;
        theta_x = 2*pi*cor_d_x/im_w;
        theta_y = pi*cor_d_y/im_h;
        theta_y_ = pi*cor_d_y_/im_h;
        r = abs(cot(theta_y)); %*ch
        cor_d_X = r*cos(-theta_x+pi/2);% *ch
        cor_d_Z = r*sin(-theta_x+pi/2);% *ch
        cor_d_Y_ = r*tan(theta_y_);% +ch
        % check correctness
        cor_a = [cor_a_X, 0, cor_a_Z];
        cor_b = [cor_b_X, 0, cor_b_Z]; 
        cor_c = [cor_c_X, 0, cor_c_Z];
        cor_d = [cor_d_X, 0, cor_d_Z];
        %keyboard
        box_w = (norm(cor_a - cor_b) + norm(cor_c - cor_d))/2;% width
        box_l = (norm(cor_b - cor_c) + norm(cor_a - cor_d))/2;% length
        % height
        box_h = mean([cor_a_Y_, cor_b_Y_, cor_c_Y_, cor_d_Y_]);
        % rotation
        cor_cen = mean([cor_a;cor_b;cor_c;cor_d]);
        cor_mid = mean([cor_a;cor_b]); 
        cor_r = cor_mid - cor_cen; cor_r = cor_r/norm(cor_r);
        cor_rc = [0,0,-1];
        box_r = acos(dot(cor_r, cor_rc));
        % normalize the box
        box_n = [box_l, box_w, box_h, cor_cen(1), cor_cen(3), box_r];
        box_n_all(rot+1,1:6) = box_n;

        %keyboard
        % flip
        cor_id(:,1) = im_w - cor_id(:,1);
        [~,corid] = sort(cor_id(:,1));
        cor_id = cor_id(corid,:);
        cor_a_x = cor_id(2,1) - im_w/2;
        cor_a_y = (im_h - cor_id(2,2))-im_h/2;
        cor_a_y_ = (im_h - cor_id(1,2))-im_h/2;
        theta_x = 2*pi*cor_a_x/im_w;
        theta_y = pi*cor_a_y/im_h;
        theta_y_ = pi*cor_a_y_/im_h;
        r = abs(cot(theta_y)); %*ch
        cor_a_X = r*cos(-theta_x+pi/2);% *ch
        cor_a_Z = r*sin(-theta_x+pi/2);% *ch 
        cor_a_Y_ = r*tan(theta_y_);% +ch

        cor_b_x = cor_id(4,1) - im_w/2;
        cor_b_y = (im_h - cor_id(4,2))-im_h/2;
        cor_b_y_ = (im_h - cor_id(3,2))-im_h/2;
        theta_x = 2*pi*cor_b_x/im_w;
        theta_y = pi*cor_b_y/im_h;
        theta_y_ = pi*cor_b_y_/im_h;
        r = abs(cot(theta_y)); %*ch
        cor_b_X = r*cos(-theta_x+pi/2);% *ch
        cor_b_Z = r*sin(-theta_x+pi/2);% *ch 
        cor_b_Y_ = r*tan(theta_y_);% +ch

        cor_c_x = cor_id(6,1) - im_w/2;
        cor_c_y = (im_h - cor_id(6,2))-im_h/2;
        cor_c_y_ = (im_h - cor_id(5,2))-im_h/2;
        theta_x = 2*pi*cor_c_x/im_w;
        theta_y = pi*cor_c_y/im_h;
        theta_y_ = pi*cor_c_y_/im_h;
        r = abs(cot(theta_y)); %*ch
        cor_c_X = r*cos(-theta_x+pi/2);% *ch
        cor_c_Z = r*sin(-theta_x+pi/2);% *ch 
        cor_c_Y_ = r*tan(theta_y_);% +ch

        cor_d_x = cor_id(8,1) - im_w/2;
        cor_d_y = (im_h - cor_id(8,2))-im_h/2;
        cor_d_y_ = (im_h - cor_id(7,2))-im_h/2;
        theta_x = 2*pi*cor_d_x/im_w;
        theta_y = pi*cor_d_y/im_h;
        theta_y_ = pi*cor_d_y_/im_h;
        r = abs(cot(theta_y)); %*ch
        cor_d_X = r*cos(-theta_x+pi/2);% *ch
        cor_d_Z = r*sin(-theta_x+pi/2);% *ch
        cor_d_Y_ = r*tan(theta_y_);% +ch
        % check correctness
        cor_a = [cor_a_X, 0, cor_a_Z];
        cor_b = [cor_b_X, 0, cor_b_Z];
        cor_c = [cor_c_X, 0, cor_c_Z];
        cor_d = [cor_d_X, 0, cor_d_Z];
        %keyboard
        box_w = (norm(cor_a - cor_b) + norm(cor_c - cor_d))/2;% width
        box_l = (norm(cor_b - cor_c) + norm(cor_a - cor_d))/2;% length
        % height
        box_h = mean([cor_a_Y_, cor_b_Y_, cor_c_Y_, cor_d_Y_]);
        % rotation
        cor_cen = mean([cor_a;cor_b;cor_c;cor_d]);
        cor_mid = mean([cor_a;cor_b]);
        cor_r = cor_mid - cor_cen; cor_r = cor_r/norm(cor_r);
        cor_rc = [0,0,-1];
        box_r = acos(dot(cor_r, cor_rc));
    
        % normalize the box
        box_n = [box_l, box_w, box_h, cor_cen(1), cor_cen(3), box_r];
        box_n_all(rot+1,7:12) = box_n;
        %keyboard
    end

    % corner map
    if 1
        kpmap = zeros(im_h,im_w);
        kpcol = zeros(1,im_w);
        for j = 1:2:size(uvcoord,1)
            kpmap_tmp = zeros(im_h,im_w);
            line = [uvcoord(j,:);uvcoord(j+1,:)];% set medium
            vert = round(mean(line(:,1)));
            line(:,1) = [vert;vert];
            [~,uvid] = sort(line(:,2));
            line = line(uvid,:);
            kpmap_tmp(line(1,2):line(2,2),vert) = 1;
            kpcol(:,vert) = 1;
            kpmap_tmp = imdilate(kpmap_tmp, se1);
            kpmap_tmp = kpmap_tmp/max(kpmap_tmp(:));
            kpmap =kpmap + kpmap_tmp;
        end
        kpmap_msk = kpmap;
        kpmap = imfilter(kpmap, h);
        kpmap = kpmap/max(kpmap(:));
        kpmap = min(kpmap,1);
        kpmap_c = kpmap;
    end

    % edge map
    if 1
        % wall
        for j = 1:2:size(cor_id,1)
            cor_id_t = cor_id(j:j+1,:);
            [ uv ] = coords2uv( cor_id_t, im_w, im_h );
            [ xyz ] = uv2xyzN( uv);
            [ lines ] = lineFromTwoPoint( xyz(1:2:end,:), xyz(2:2:end,:) );
            [ panoEdgeC ] = paintParameterLine_my(lines, im_w, im_h, panoImg);
            panoEdgeC = panoEdgeC(:,:,1);
            panoEdgeC = imdilate(panoEdgeC, [se1, se2]);
            panoEdgeC = imfilter(panoEdgeC, h);
            panoEdgeC = panoEdgeC/max(panoEdgeC(:));
            kpmap(:,:,1) = max(kpmap(:,:,1), panoEdgeC);
        end
    
        % ceiling
        cor_all = [cor_id(1,:);cor_id(3,:);cor_id(3,:);cor_id(5,:);
               cor_id(5,:);cor_id(7,:);cor_id(7,:);cor_id(1,:)];
        for j = 1:2:size(cor_all, 1)
            cor_id_t = cor_all(j:j+1,:);
            [ uv ] = coords2uv( cor_id_t, im_w, im_h );
            [ xyz ] = uv2xyzN( uv);
            [ lines ] = lineFromTwoPoint( xyz(1:2:end,:), xyz(2:2:end,:) );
            [ panoEdgeC ] = paintParameterLine_my(lines, im_w, im_h, panoImg);
            panoEdgeC = panoEdgeC(:,:,1);
            panoEdgeC = imdilate(panoEdgeC, [se1, se2]);
            panoEdgeC = imfilter(panoEdgeC, h);
            panoEdgeC = panoEdgeC/max(panoEdgeC(:));
            kpmap(:,:,2) = max(kpmap(:,:,2), panoEdgeC);
        end
    
        cor_all = [cor_id(2,:);cor_id(4,:);cor_id(4,:);cor_id(6,:);
                   cor_id(6,:);cor_id(8,:);cor_id(8,:);cor_id(2,:)];
        % floor       
        for j = 1:2:size(cor_all, 1)
            cor_id_t = cor_all(j:j+1,:);
            [ uv ] = coords2uv( cor_id_t, im_w, im_h );
            [ xyz ] = uv2xyzN( uv);
            [ lines ] = lineFromTwoPoint( xyz(1:2:end,:), xyz(2:2:end,:) );
            [ panoEdgeC ] = paintParameterLine_my(lines, im_w, im_h, panoImg);
            panoEdgeC = panoEdgeC(:,:,1);
            panoEdgeC = imdilate(panoEdgeC, [se1, se2]);
            panoEdgeC = imfilter(panoEdgeC, h);
            panoEdgeC = panoEdgeC/max(panoEdgeC(:));
            kpmap(:,:,3) = max(kpmap(:,:,3), panoEdgeC);
        end
    
        kpmap = min(kpmap,1);
    
        %keyboard
        % visual
        if visual
            pano_vis = panoImg;
            pano_vis(:,:,1) = pano_vis(:,:,1).*kpmap_msk;
            ha = tight_subplot(3, 1);
            axes(ha(1));imagesc(panoImg);axis image;
            axes(ha(2));imagesc(pano_vis);axis image; 
            axes(ha(3));imagesc(kpmap);axis image;
            keyboard
        end
    end
    
    save(['../data/pano_edge_tr/box/' im_name '.mat' ], 'box_n_all');
    save(['../data/pano_edge_tr/vp/' im_name '.mat' ], 'vp');
    imwrite(kpmap, ['../data/pano_edge_tr/lay/' im_name '.png' ],'png');
    imwrite(kpmap_c, ['../data/pano_edge_tr/cor/' im_name '.png' ],'png');
    imwrite(panoImg, ['../data/pano_edge_tr/img/' im_name '.png' ],'png');
end
