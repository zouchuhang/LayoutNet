
% path
addpath(genpath('./panoContext_code/'));
addpath(genpath('./minFunc_2012/'));

% LBFGS options
options = [];
options.display = 'none';
options.MaxIter = 100;
options.Method = 'lbfgs';
options.LS_init = 2;

% test image name
fid = fopen('../gt/panoContext_test.txt');
tline = fgetl(fid);

visual = 1;

% visualization
if visual
    figure(1);
    ha = tight_subplot(2,2);
end

for i = 1:53

    disp(tline);

    % read LaoutNet results
    im = imread(['../result/res_panofull_ts_box_joint/img/' num2str(i) '.png']);
    edg = imread(['../result/res_panofull_ts_box_joint/edg/' num2str(i) '.png']);
    corn = imread(['../result/res_panofull_ts_box_joint/cor/' num2str(i) '.png']);
    
    edg2 = edg;
    edg = edg(:,:,1);
    % max
    edg_m = max(edg); edg_m = double(edg_m);
    cor_m = max(corn); cor_m = double(cor_m);

    [im_h, im_w, ~] = size(im);
    
    % find peak response
    [cor_id, pks, pk_loc] = getIniCor(cor_m, corn, edg_m, im_h);
    
    % sampling
    disp('sampling ...')
    tic
    cor_fn = samplingPanoBox(cor_id, corn, edg, edg2, im_h, im_w, options);
    toc
    cor_id = cor_fn;

    tline = fgetl(fid);

    if visual
        % visual 
        axes(ha(1));cla;
        imagesc(im); axis image; hold on;
        % rescale
        edg_m = (edg_m - min(edg_m))/(max(edg_m) - min(edg_m)) * 100;
        edg_m = edg_m + 256-100;
        plot(1:im_w, edg_m, 'r', 'LineWidth',2);
        for j = 1:min(4, numel(pks))
            plot([pk_loc(j) pk_loc(j)], [0, im_h], '--b', 'LineWidth',1);
        end
        title('Wall-wall Boundaries');
        hold off; 

        cor_all = [cor_id;
               cor_id(1,:);cor_id(3,:);cor_id(3,:);cor_id(5,:);
               cor_id(5,:);cor_id(7,:);cor_id(7,:);cor_id(1,:);
               cor_id(2,:);cor_id(4,:);cor_id(4,:);cor_id(6,:);
               cor_id(6,:);cor_id(8,:);cor_id(8,:);cor_id(2,:)];
        [ uv ] = coords2uv( cor_all, im_w, im_h );
        [ xyz ] = uv2xyzN( uv);
        [ lines ] = lineFromTwoPoint( xyz(1:2:end,:), xyz(2:2:end,:) );
        [ panoEdgeC ] = paintParameterLine2(lines, im_w, im_h, im, 1 );
        axes(ha(2));
        imagesc(panoEdgeC); title('Optimized Layout'); axis image; hold on;

        % corners
        for j = 1:size(cor_id,1)
            plot(cor_id(j,1),cor_id(j,2), '+b');
        end
        hold off; 
        axes(ha(3));
        imagesc(edg); title('Edge Map'); axis image;
        hold off; 
        axes(ha(4));
        imagesc(corn); title('Corner Map'); axis image;

        keyboard
    end
end
