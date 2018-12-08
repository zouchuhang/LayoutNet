
addpath(genpath('./LSUN/toolkit'))
data_path = './LSUN/';
val_set = load([data_path '/test.mat']);
% test list
nl = val_set.test;
name_list = {nl.image};

% match test images in .t7 file with the gt name for evaluation
fid = fopen('../data/lsun_ts.txt');
tline = fgetl(fid);

figure(1);ha = tight_subplot(2,2);

RoomLayoutTypes;
type(8).cornermap = [1 3 7 5];

im_h = 512;
im_w = 512;

for i = 1:1000
    
    
    disp(i);

    im = imread(['../data/lsun_ts/img/' tline]);

    edg = imread(['../result/res_lsun_ts_512_joint/edg/' num2str(i) '.png']);
    corn = load(['../result/res_lsun_ts_512_joint/cor_mat/' num2str(i) '.mat']); corn = corn.x;
    corn = permute(corn,[2,3,1]);
    corn_f = load(['../result/res_lsun_ts_512_joint/cor_mat_flip/' num2str(i) '.mat']); corn_f = corn_f.x;
    corn_f = permute(corn_f,[2,3,1]);

    % find room type
    if 1
        r_t = load(['../result/res_lsun_ts_512_joint/type/' num2str(i) '.mat']); r_t = r_t.x;
        r_t = mean(r_t);
        [~,RecordId] = max(r_t);
    end

    if 0
    % assume known type
        RecordId = find([type.typeid] == nl(idx).type);
    end

    room_t = type(RecordId);
    im_ori = imread([data_path 'image/images/' tline(1:end-4) '.jpg']);
    %im_res = nl(idx).resolution;
    im_res(1) = size(im_ori,1);
    im_res(2) = size(im_ori,2);


    % flip
    if room_t.typeid == 0
        corn_t = corn_f(:,:,1); corn_f(:,:,1) = corn_f(:,:,7); corn_f(:,:,7) = corn_t;
        corn_t = corn_f(:,:,3); corn_f(:,:,3) = corn_f(:,:,5); corn_f(:,:,5) = corn_t;
        corn_t = corn_f(:,:,2); corn_f(:,:,2) = corn_f(:,:,8); corn_f(:,:,8) = corn_t;
        corn_t = corn_f(:,:,4); corn_f(:,:,4) = corn_f(:,:,6); corn_f(:,:,6) = corn_t;
        corn(:,:,3) = max(0,corn(:,:,3) - corn(:,:,4));
        corn(:,:,5) = max(0,corn(:,:,5) - corn(:,:,6));
        corn(:,:,1) = max(0,corn(:,:,1) - corn(:,:,2));
        corn(:,:,7) = max(0,corn(:,:,7) - corn(:,:,8));
        corn_f(:,:,3) = max(0,corn_f(:,:,3) - corn_f(:,:,4));
        corn_f(:,:,5) = max(0,corn_f(:,:,5) - corn_f(:,:,6));
        corn_f(:,:,1) = max(0,corn_f(:,:,1) - corn_f(:,:,2));
        corn_f(:,:,7) = max(0,corn_f(:,:,7) - corn_f(:,:,8));
        corn(:,400:end,1) = 0; corn(:,400:end,3) = 0;
        corn(:,1:112,7) = 0; corn(:,1:112,7) = 0;
        
        %keyboard
        a = corn(:,:,5); b = corn(:,:,3); corn(:,:,5)=max(corn(:,:,5)-b, 0); corn(:,:,3)=max(corn(:,:,3)-a,0);
        a = corn(:,:,1); b = corn(:,:,7); corn(:,:,1)=max(corn(:,:,1)-b,0); corn(:,:,7)=max(corn(:,:,7)-a,0);
        a = corn_f(:,:,5); b = corn_f(:,:,3); corn_f(:,:,5)=max(corn_f(:,:,5)-b,0); corn_f(:,:,3)=max(corn_f(:,:,3)-a,0);
        a = corn_f(:,:,1); b = corn_f(:,:,7); corn_f(:,:,1)=max(corn_f(:,:,1)-b,0); corn_f(:,:,7)=max(corn_f(:,:,7)-a,0);
    end
    if room_t.typeid == 1
        corn_t = corn_f(:,:,1); corn_f(:,:,1) = corn_f(:,:,7); corn_f(:,:,7) = corn_t;
        corn_t = corn_f(:,:,3); corn_f(:,:,3) = corn_f(:,:,5); corn_f(:,:,5) = corn_t;
        corn_t = corn_f(:,:,4); corn_f(:,:,4) = corn_f(:,:,6); corn_f(:,:,6) = corn_t;
        corn(:,:,3) = max(0,corn(:,:,3) - corn(:,:,4));
        corn(:,:,5) = max(0,corn(:,:,5) - corn(:,:,6));
        corn_f(:,:,3) = max(0,corn_f(:,:,3) - corn_f(:,:,4));
        corn_f(:,:,5) = max(0,corn_f(:,:,5) - corn_f(:,:,6));
        
        a = corn(:,:,5); b = corn(:,:,3); corn(:,:,5)=max(corn(:,:,5)-b, 0); corn(:,:,3)=max(corn(:,:,3)-a,0);
        a = corn(:,:,1); b = corn(:,:,7); corn(:,:,1)=max(corn(:,:,1)-b,0); corn(:,:,7)=max(corn(:,:,7)-a,0);
        a = corn_f(:,:,5); b = corn_f(:,:,3); corn_f(:,:,5)=max(corn_f(:,:,5)-b,0); corn_f(:,:,3)=max(corn_f(:,:,3)-a,0);
        a = corn_f(:,:,1); b = corn_f(:,:,7); corn_f(:,:,1)=max(corn_f(:,:,1)-b,0); corn_f(:,:,7)=max(corn_f(:,:,7)-a,0);
    end
    if room_t.typeid == 2
        corn_t = corn_f(:,:,2); corn_f(:,:,2) = corn_f(:,:,8); corn_f(:,:,8) = corn_t;
        corn_t = corn_f(:,:,1); corn_f(:,:,1) = corn_f(:,:,7); corn_f(:,:,7) = corn_t;
        corn_t = corn_f(:,:,3); corn_f(:,:,3) = corn_f(:,:,5); corn_f(:,:,5) = corn_t;
        corn(:,:,1) = max(0,corn(:,:,1) - corn(:,:,2));
        corn(:,:,7) = max(0,corn(:,:,7) - corn(:,:,8));
        corn_f(:,:,1) = max(0,corn_f(:,:,1) - corn_f(:,:,2));
        corn_f(:,:,7) = max(0,corn_f(:,:,7) - corn_f(:,:,8));
        
        a = corn(:,:,5); b = corn(:,:,3); corn(:,:,5)=max(corn(:,:,5)-b, 0); corn(:,:,3)=max(corn(:,:,3)-a,0);
        a = corn(:,:,1); b = corn(:,:,7); corn(:,:,1)=max(corn(:,:,1)-b,0); corn(:,:,7)=max(corn(:,:,7)-a,0);
        a = corn_f(:,:,5); b = corn_f(:,:,3); corn_f(:,:,5)=max(corn_f(:,:,5)-b,0); corn_f(:,:,3)=max(corn_f(:,:,3)-a,0);
        a = corn_f(:,:,1); b = corn_f(:,:,7); corn_f(:,:,1)=max(corn_f(:,:,1)-b,0); corn_f(:,:,7)=max(corn_f(:,:,7)-a,0);
    end
    if room_t.typeid == 3
        corn_t = corn_f(:,:,1); corn_f(:,:,1) = corn_f(:,:,8);corn_f(:,:,8) = corn_t;
        corn(:,:,7) = max(0,corn(:,:,7) - corn(:,:,1));
        corn_f(:,:,7) = max(0,corn_f(:,:,7) - corn_f(:,:,1));
        corn(:,:,7) = max(0,corn(:,:,7) - corn(:,:,8));
        corn_f(:,:,7) = max(0,corn_f(:,:,7) - corn_f(:,:,8));
        corn(:,:,5) = max(0,corn(:,:,5) - corn(:,:,3));
        corn_f(:,:,5) = max(0,corn_f(:,:,5) - corn_f(:,:,3));
    end
    if room_t.typeid == 4
        corn_t = corn_f(:,:,3); corn_f(:,:,3) = corn_f(:,:,6); corn_f(:,:,6) = corn_t;
        corn(:,:,5) = max(0,corn(:,:,5) - corn(:,:,3));
        corn_f(:,:,5) = max(0,corn_f(:,:,5) - corn_f(:,:,3));
        corn(:,:,5) = max(0,corn(:,:,5) - corn(:,:,6));
        corn_f(:,:,5) = max(0,corn_f(:,:,5) - corn_f(:,:,6));
        corn(:,:,7) = max(0,corn(:,:,7) - corn(:,:,1));
        corn_f(:,:,7) = max(0,corn_f(:,:,7) - corn_f(:,:,1));
    end
    if room_t.typeid == 5
        corn_t = corn_f(:,:,1); corn_f(:,:,1) = corn_f(:,:,8);corn_f(:,:,8) = corn_t;
        corn_t = corn_f(:,:,3); corn_f(:,:,3) = corn_f(:,:,6);corn_f(:,:,6) = corn_t;
        corn(:,:,7) = max(0,corn(:,:,7) - corn(:,:,1));
        corn_f(:,:,7) = max(0,corn_f(:,:,7) - corn_f(:,:,1));
        corn(:,:,5) = max(0,corn(:,:,5) - corn(:,:,3));
        corn_f(:,:,5) = max(0,corn_f(:,:,5) - corn_f(:,:,3));
        corn(:,:,7) = max(0,corn(:,:,7) - corn(:,:,8));
        corn_f(:,:,7) = max(0,corn_f(:,:,7) - corn_f(:,:,8));
        corn(:,:,5) = max(0,corn(:,:,5) - corn(:,:,6));
        corn_f(:,:,5) = max(0,corn_f(:,:,5) - corn_f(:,:,6));
    end
    if room_t.typeid == 6
        corn_t = corn_f(:,:,1); corn_f(:,:,1) = corn_f(:,:,7);corn_f(:,:,7) = corn_t;
        corn_t = corn_f(:,:,3); corn_f(:,:,3) = corn_f(:,:,5);corn_f(:,:,5) = corn_t;
    end
    if room_t.typeid == 7
        keyboard
    end
    if room_t.typeid == 8
        keyboard
    end
    if room_t.typeid == 9
        corn_t = corn_f(:,:,3); corn_f(:,:,3) = corn_f(:,:,5);corn_f(:,:,5) = corn_t;
    end

    % find corner
    point = [];
    for j = 1:numel(room_t.cornermap)
        mp = corn(:,:,room_t.cornermap(j)) + fliplr(corn_f(:,:,room_t.cornermap(j)));
        mp(:,1) = 0;mp(:,im_w) = 0; mp(1,:) = 0;mp(im_h,:) = 0;
        %mp = mp(1:size(mp,1)-1, 1:size(mp,2)-1);
        mp_msk = zeros(size(mp));
        if room_t.typeid == 1
            if room_t.cornermap(j) == 7 || room_t.cornermap(j) == 1
                mp_msk = edg(:,:,2) >255*0.1;
            end
            if room_t.cornermap(j) == 3
                mp_msk = edg(:,:,2) >255*0.1;
                mp_t = (corn(:,:,1) + fliplr(corn_f(:,:,1)))/2;
                mp_t(:,1) = 0;mp_t(:,im_w) = 0; mp_t(1,:) = 0;mp_t(im_h,:) = 0;
                mp_t = mp_t.*mp_msk;
                [~,pt] = max(mp_t(:));
                [~, pt_x] = ind2sub([im_h im_w], pt);
                mp_msk(:, 1:max(pt_x - 50,1)) = 0;
                mp_msk(:, min(pt_x + 50,im_w):end) = 0;
            end
            if room_t.cornermap(j) == 5
                mp_msk = edg(:,:,2) >255*0.1;
                mp_t = (corn(:,:,7) + fliplr(corn_f(:,:,7)))/2;
                mp_t(:,1) = 0;mp_t(:,im_w) = 0; mp_t(1,:) = 0;mp_t(im_h,:) = 0;
                mp_t = mp_t.*mp_msk;
                [~,pt] = max(mp_t(:));
                [~, pt_x] = ind2sub([im_h im_w], pt);
                mp_msk(:, 1:max(pt_x - 50,1)) = 0;
                mp_msk(:, min(pt_x + 50,im_w):end) = 0;
            end
            if room_t.cornermap(j) == 4 || room_t.cornermap(j) == 6
                mp_msk = edg(:,:,3) >255*0.1;
            end
        end
        if room_t.typeid == 2
            keyboard
            if room_t.cornermap(j) == 2 || room_t.cornermap(j) == 8
                mp_msk = edg(:,:,1) >255*0.1;
            end
            if room_t.cornermap(j) == 5 || room_t.cornermap(j) == 3 || room_t.cornermap(j) == 1 || room_t.cornermap(j) == 7
                mp_msk = edg(:,:,2) >255*0.1;
            end
        end
        if room_t.typeid == 3
            if room_t.cornermap(j) == 1 || room_t.cornermap(j) == 8
                mp_msk = edg(:,:,1) >255*0.1;
            end
            if room_t.cornermap(j) == 5 || room_t.cornermap(j) == 7
                mp_msk = edg(:,:,2) >255*0.1;
            end
        end
        if room_t.typeid == 4
            if room_t.cornermap(j) == 5
                mp_msk = edg(:,:,2) >255*0.1;
            end
            if room_t.cornermap(j) == 7
                mp_msk = edg(:,:,2) >255*0.1;
                mp_t = (corn(:,:,5) + fliplr(corn_f(:,:,5)))/2;
                mp_t(:,1) = 0;mp_t(:,im_w) = 0; mp_t(1,:) = 0;mp_t(im_h,:) = 0;
                mp_t = mp_t.*mp_msk;
                [~,pt] = max(mp_t(:));
                [~, pt_x] = ind2sub([im_h im_w], pt);
                mp_msk(:, 1:max(pt_x - 50,1)) = 0;
                mp_msk(:, min(pt_x + 50,im_w):end) = 0;
            end 
            if room_t.cornermap(j) == 3 || room_t.cornermap(j) == 6
                mp_msk = edg(:,:,3) >255*0.1;
            end
        end
        if room_t.typeid == 5
            if room_t.cornermap(j) == 7 %|| room_t.cornermap(j) == 5
                mp_msk = edg(:,:,2) >255*0.1;
            end
            if room_t.cornermap(j) == 5
                mp_msk = edg(:,:,2) >255*0.1;
                mp_t = (corn(:,:,7) + fliplr(corn_f(:,:,7)))/2;
                mp_t(:,1) = 0;mp_t(:,im_w) = 0; mp_t(1,:) = 0;mp_t(im_h,:) = 0;
                mp_t = mp_t.*mp_msk;
                [~,pt] = max(mp_t(:));
                [~, pt_x] = ind2sub([im_h im_w], pt);
                mp_msk(:, 1:max(pt_x - 50,1)) = 0;
                mp_msk(:, min(pt_x + 50,im_w):end) = 0;
            end
            if room_t.cornermap(j) == 1 || room_t.cornermap(j) == 8
                mp_msk = edg(:,:,1) >255*0.1;
            end
            if room_t.cornermap(j) == 3 || room_t.cornermap(j) == 6
                mp_msk = edg(:,:,3) >255*0.1;
            end
        end
        if room_t.typeid == 6
            %keyboard
            if room_t.cornermap(j) == 1 || room_t.cornermap(j) == 7
                mp_msk = edg(:,:,1) >255*0.1;
            end
            if room_t.cornermap(j) == 3 || room_t.cornermap(j) == 5
                mp_msk = edg(:,:,3) >255*0.1;
            end
        end
        if room_t.typeid == 9
            if room_t.cornermap(j) == 5 || room_t.cornermap(j) == 3
                mp_msk = edg(:,:,3) >255*0.1;
            end
        end
         if room_t.typeid == 10
            if room_t.cornermap(j) == 5 || room_t.cornermap(j) == 7
                mp_msk = edg(:,:,2) >255*0.1;
            end
        end
        if room_t.typeid == 0
            if room_t.cornermap(j) == 5 || room_t.cornermap(j) == 3
                mp_msk = edg(:,:,2) >255*0.1;
            end

            if room_t.cornermap(j) == 1 %|| room_t.cornermap(j) == 7
                mp_msk = edg(:,:,2) >255*0.1;
                mp_t = (corn(:,:,3) + fliplr(corn_f(:,:,3)))/2;
                mp_t(:,1) = 0;mp_t(:,im_w) = 0; mp_t(1,:) = 0;mp_t(im_h,:) = 0;
                mp_t = mp_t.*mp_msk;
                [~,pt] = max(mp_t(:));
                [~, pt_x] = ind2sub([im_h im_w], pt);
                mp_msk(:, 1:max(pt_x - 50,1)) = 0;
                mp_msk(:, min(pt_x + 50,im_w):end) = 0;
            end
            if room_t.cornermap(j) == 7
                mp_msk = edg(:,:,2) >255*0.1;
                mp_t = (corn(:,:,5) + fliplr(corn_f(:,:,5)))/2;
                mp_t(:,1) = 0;mp_t(:,im_w) = 0; mp_t(1,:) = 0;mp_t(im_h,:) = 0;
                mp_t = mp_t.*mp_msk;
                [~,pt] = max(mp_t(:));
                [~, pt_x] = ind2sub([im_h im_w], pt);
                mp_msk(:, 1:max(pt_x - 50,1)) = 0;
                mp_msk(:, min(pt_x + 50,im_w):end) = 0;
            end
            if room_t.cornermap(j) == 8 || room_t.cornermap(j) == 2
                mp_msk = edg(:,:,1) >255*0.1;
            end
            if room_t.cornermap(j) == 4 || room_t.cornermap(j) == 6
                mp_msk = edg(:,:,3) >255*0.1;
            end
        end
        mp = mp.*mp_msk;
        
        [~,pt] = max(mp(:)/2);
        [pt_x, pt_y] = ind2sub([im_h im_w], pt);
        point = [point;pt_x pt_y];
        %keyboard
    end

    point_res = [point(:,2) point(:,1)];
    cor_res = im_res./[im_w im_h];
    point_ref_res = bsxfun(@times, point_res, [cor_res(2) cor_res(1)]);
    point_ref_res(point_ref_res(:,1) > im_res(2),1) = im_res(2);
    point_ref_res(point_ref_res(:,2) > im_res(1),2) = im_res(1);
    point_ref_res(point_ref_res<1) = 1;
    %keyboard
    % refine point
    point = point_ref_res;
    point_ref = point;
    P = [0 0; 0 im_res(1)+0.01; im_res(2)+0.01 im_res(1)+0.01; im_res(2)+0.01 0]; 
    P = P'; P = [P P(:,1)];
    if room_t.typeid == 0
        line_1 = polyfit([point(1,1) point(2,1)], [point(1,2) point(2,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(1,1); point(1,2)];
        s1(:,2) = [-100; -100*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(2,:) = X';
        line_1 = polyfit([point(3,1) point(4,1)], [point(3,2) point(4,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(3,1); point(3,2)];
        s1(:,2) = [-100; -100*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(4,:) = X';
        line_1 = polyfit([point(5,1) point(6,1)], [point(5,2) point(6,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(5,1); point(5,2)];
        s1(:,2) = [100000; 100000*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(6,:) = X';
        line_1 = polyfit([point(7,1) point(8,1)], [point(7,2) point(8,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(7,1); point(7,2)];
        s1(:,2) = [100000; 100000*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(8,:) = X';
    end
    if room_t.typeid == 1
        line_1 = polyfit([point(1,1) point(2,1)], [point(1,2) point(2,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(1,1); point(1,2)];
        s1(:,2) = [(-1-line_1(2))/line_1(1);-1];
        X = seg2poly(s1, P); point_ref(2,:) = X';
        line_1 = polyfit([point(1,1) point(3,1)], [point(1,2) point(3,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(1,1); point(1,2)];
        s1(:,2) = [-100; -100*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(3,:) = X';
        if point(4,1) == point(5,1)
            point(5,1) = point(5,1) + 0.01;
        end
        line_1 = polyfit([point(4,1) point(5,1)], [point(4,2) point(5,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(4,1); point(4,2)];
        s1(:,2) = [(-1-line_1(2))/line_1(1);-1];
        X = seg2poly(s1, P); point_ref(5,:) = X';
        line_1 = polyfit([point(4,1) point(6,1)], [point(4,2) point(6,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(4,1); point(4,2)];
        s1(:,2) = [100000; 100000*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(6,:) = X';
    end
    if room_t.typeid == 2
        keyboard
    end
    if room_t.typeid == 3
        line_1 = polyfit([point(1,1) point(2,1)], [point(1,2) point(2,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(1,1); point(1,2)];
        s1(:,2) = [-100; -100*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(2,:) = X';
        line_1 = polyfit([point(1,1) point(4,1)], [point(1,2) point(4,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(1,1); point(1,2)];
        s1(:,2) = [100000; 100000*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(4,:) = X';
        line_1 = polyfit([point(1,1) point(3,1)], [point(1,2) point(3,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(1,1); point(1,2)];
        s1(:,2) = [(10000-line_1(2))/line_1(1);10000];
        X = seg2poly(s1, P); point_ref(3,:) = X';
    end
    if room_t.typeid == 4
        line_1 = polyfit([point(1,1) point(2,1)], [point(1,2) point(2,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(1,1); point(1,2)];
        s1(:,2) = [-100; -100*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(2,:) = X';
        if point(1,1) == point(3,1)
            point(1,1) = point(1,1) + 0.01;
        end
        line_1 = polyfit([point(1,1) point(3,1)], [point(1,2) point(3,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(1,1); point(1,2)];
        s1(:,2) = [(-1-line_1(2))/line_1(1);-1];
        X = seg2poly(s1, P); point_ref(3,:) = X';
        line_1 = polyfit([point(1,1) point(4,1)], [point(1,2) point(4,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(1,1); point(1,2)];
        s1(:,2) = [10000; 10000*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(4,:) = X';
    end
    if room_t.typeid == 5
        line_1 = polyfit([point(1,1) point(2,1)], [point(1,2) point(2,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(1,1); point(1,2)];
        s1(:,2) = [-100; -100*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(2,:) = X';
        line_1 = polyfit([point(1,1) point(3,1)], [point(1,2) point(3,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(1,1); point(1,2)];
        s1(:,2) = [10000; 10000*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(3,:) = X';
        line_1 = polyfit([point(4,1) point(5,1)], [point(4,2) point(5,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(4,1); point(4,2)];
        s1(:,2) = [-100; -100*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(5,:) = X';
        line_1 = polyfit([point(4,1) point(6,1)], [point(4,2) point(6,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(4,1); point(4,2)];
        s1(:,2) = [10000; 10000*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(6,:) = X';
    end
    if room_t.typeid == 6
        line_1 = polyfit([point(1,1) point(2,1)], [point(1,2) point(2,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(1,1); point(1,2)];
        s1(:,2) = [-100; -100*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(1,:) = X';
        s1 = zeros(2,2); s1(:,1) = [point(2,1); point(2,2)];
        s1(:,2) = [10000; 10000*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(2,:) = X';
        line_1 = polyfit([point(3,1) point(4,1)], [point(3,2) point(4,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(3,1); point(3,2)];
        s1(:,2) = [-100; -100*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(3,:) = X';
        s1 = zeros(2,2); s1(:,1) = [point(4,1); point(4,2)];
        s1(:,2) = [10000; 10000*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(4,:) = X';
    end
    if room_t.typeid == 7
        keyboard
    end
    if room_t.typeid == 8
        keyboard
    end
    if room_t.typeid == 9
        line_1 = polyfit([point(1,1) point(2,1)], [point(1,2) point(2,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(1,1); point(1,2)];
        s1(:,2) = [-100; -100*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(1,:) = X';
        s1 = zeros(2,2); s1(:,1) = [point(2,1); point(2,2)];
        s1(:,2) = [10000; 10000*line_1(1) + line_1(2)];
        X = seg2poly(s1, P); point_ref(2,:) = X';
    end
    if room_t.typeid == 10
        line_1 = polyfit([point(1,1) point(2,1)], [point(1,2) point(2,2)], 1);
        s1 = zeros(2,2); s1(:,1) = [point(1,1); point(1,2)];
        s1(:,2) = [(-1-line_1(2))/line_1(1);-1];
        X = seg2poly(s1, P); point_ref(1,:) = X';
        s1 = zeros(2,2); s1(:,1) = [point(2,1); point(2,2)];
        s1(:,2) = [(10000-line_1(2))/line_1(1);10000];
        X = seg2poly(s1, P); point_ref(2,:) = X';
    end

    % get segment
    data.type = room_t.typeid; data.point = point_ref; data.resolution = im_res;
    [ seg ] = getSegmentation( data );
    %[gt_seg] = getSegmentation( nl(idx) );
    data.layout = seg;

    % acc
    %PtError = cornerError(point_ref, nl(idx).point, im_res);
    %allPxError = 1 - pixelwiseAccuracy(seg, gt_seg, im_res);
    %acc_co_all(i) = PtError;
    %lay_acc_all(i) = allPxError;
    % vis
    axes(ha(1));imagesc(im_ori); axis image; %title(PtError);
    axes(ha(2));imagesc(data.layout); axis image; %title(allPxError);
    axes(ha(3));imagesc(imresize(edg, data.resolution)); axis image;
    axes(ha(4));imagesc(imresize(sum(corn, 3), data.resolution)); axis image;

    result = data;

    save(['../result/res_lsun_ts_512_joint/mat_ts_v73/' tline(1:end-4) '.mat'], 'result', '-v7.3');
    %keyboard

    tline = fgetl(fid);
end
