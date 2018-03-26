function [cor_ini, score_all] = sample_opt_pano(cor_id, im_w, im_h, corn, edg, edg2, options)

	score_all = 0;
        [wall_d, x, f] = pano_line_solver(cor_id, im_w, options); % bos-shape case
        cor_id_t = cor_id;
        d_m = 0;
        % initialization
        for j = 1:2:size(cor_id_t,1)
            theta_y = pi*cor_id_t(j+1,2)/im_h-pi/2;
            theta_y_ = pi/2 - pi*cor_id_t(j,2)/im_h;
            d = abs(cot(theta_y));
            d_m = d_m + d/wall_d((j+1)/2);
        end
        d_m = d_m/4;
        cor_ini = [];
        for j = 1:2:size(cor_id_t,1)
            line_d = d_m*wall_d((j+1)/2);
            line_theta_y = acot(line_d);
            line_y = (line_theta_y + pi/2)*im_h/pi;
            cor_ini = [cor_ini;cor_id_t(j,1) line_y];
            cor_ini = [cor_ini;cor_id_t(j,1) line_y];
        end
        
        score_btn = interp2(corn,cor_ini(2:2:end,1),cor_ini(2:2:end,2));
        score_btn = sum(log(score_btn));
       
        % add floor score
        cor_all = [cor_ini(2,:);cor_ini(4,:);cor_ini(4,:);cor_ini(6,:);
               cor_ini(6,:);cor_ini(8,:);cor_ini(8,:);cor_ini(2,:)];
    	[ uv ] = coords2uv( cor_all, im_w, im_h );
    	[ xyz ] = uv2xyzN( uv);
    	[ lines ] = lineFromTwoPoint( xyz(1:2:end,:), xyz(2:2:end,:) );
    	im = zeros(im_h, im_w, 3);
    	[ panoEdgeC ] = paintParameterLine_my(lines, im_w, im_h, im);
    	score_fl = zeros(1,4);
    	score_fl(1) = max(max(edg2(:,:,3).*(panoEdgeC(:,:,1) == 1)));
    	score_fl(2) = max(max(edg2(:,:,3).*(panoEdgeC(:,:,1) == 2)));
    	score_fl(3) = max(max(edg2(:,:,3).*(panoEdgeC(:,:,1) == 3)));
    	score_fl(4) = max(max(edg2(:,:,3).*(panoEdgeC(:,:,1) == 4)));
    	score_fl = sum(log(score_fl));
    	score_btn = score_btn + score_fl;

        % sampling
        % horizontal, bottoms
        d_m_max = d_m + d_m*0.1;
        d_m_min = d_m - d_m*0.1; 
        score_best = -Inf;

        for j = d_m_min:d_m*0.02:d_m_max
            cor_opt = [];
            for k = 1:2:size(cor_id_t,1)
                line_d = j*wall_d((k+1)/2);
                line_theta_y = acot(line_d);
                line_y = (line_theta_y + pi/2)*im_h/pi;
                cor_opt = [cor_opt;cor_id_t(k,1) line_y];
            end
            % compute score
            score = interp2(corn,cor_opt(:,1),cor_opt(:,2));
            score = sum(log(score));

            cor_all = [cor_opt(1,:);cor_opt(2,:);cor_opt(2,:);cor_opt(3,:);
               cor_opt(3,:);cor_opt(4,:);cor_opt(4,:);cor_opt(1,:)];
    		[ uv ] = coords2uv( cor_all, im_w, im_h );
    		[ xyz ] = uv2xyzN( uv);
    		[ lines ] = lineFromTwoPoint( xyz(1:2:end,:), xyz(2:2:end,:) );
    		[ panoEdgeC ] = paintParameterLine_my(lines, im_w, im_h, im);
    		score_fl = zeros(1,4);
    		score_fl(1) = max(max(edg2(:,:,3).*(panoEdgeC(:,:,1) == 1)));
    		score_fl(2) = max(max(edg2(:,:,3).*(panoEdgeC(:,:,1) == 2)));
    		score_fl(3) = max(max(edg2(:,:,3).*(panoEdgeC(:,:,1) == 3)));
    		score_fl(4) = max(max(edg2(:,:,3).*(panoEdgeC(:,:,1) == 4)));
    		score_fl = sum(log(score_fl));
    		score = score + score_fl;

            if score > score_best
                score_best = score;
                cor_best = cor_opt;
                d_m_best = j;
            end
        end
        if score_best > score_btn
            cor_ini(2:2:end,:) = cor_best;
        else
           score_best = score_btn;
           d_m_best = d_m;
        end
        score_all = score_all + score_best;
        % horizontal, top
        h_m = 0;
        for j = 1:2:size(cor_id_t,1)
            theta_y_ = pi/2 - pi*cor_id_t(j,2)/im_h;
            d = d_m_best*wall_d((j+1)/2);
            h = d*tan(theta_y_);
            h_m = h_m + h;
        end
        h_m = h_m/4;
        for j = 1:2:size(cor_id_t,1)
            line_d = d_m_best*wall_d((j+1)/2);
            line_y_ = atan(h_m/line_d);
            line_y_ = (pi/2-line_y_)*im_h/pi;
            cor_ini(j,2) = line_y_;
        end
        score_tp = interp2(corn,cor_ini(1:2:end,1),cor_ini(1:2:end,2));
        score_tp = sum(log(score_tp));

        % add ceiling score
        cor_all = [cor_ini(1,:);cor_ini(3,:);cor_ini(3,:);cor_ini(5,:);
               cor_ini(5,:);cor_ini(7,:);cor_ini(7,:);cor_ini(1,:)];
    	[ uv ] = coords2uv( cor_all, im_w, im_h );
    	[ xyz ] = uv2xyzN( uv);
    	[ lines ] = lineFromTwoPoint( xyz(1:2:end,:), xyz(2:2:end,:) );
    	im = zeros(im_h, im_w, 3);
    	[ panoEdgeC ] = paintParameterLine_my(lines, im_w, im_h, im);
    	score_cl = zeros(1,4);

    	score_cl(1) = max(max(edg2(:,:,2).*(panoEdgeC(:,:,1) == 1)));
    	score_cl(2) = max(max(edg2(:,:,2).*(panoEdgeC(:,:,1) == 2)));
    	score_cl(3) = max(max(edg2(:,:,2).*(panoEdgeC(:,:,1) == 3)));
    	score_cl(4) = max(max(edg2(:,:,2).*(panoEdgeC(:,:,1) == 4)));
    	score_cl = sum(log(score_cl));
    	score_tp = score_tp + 0.5*score_cl;

        h_m_max = h_m + h_m*0.1;
        h_m_min = h_m - h_m*0.1;
        score_best = -Inf;
        for j = h_m_min:h_m*0.02:h_m_max
            cor_opt = [];
            for k = 1:2:size(cor_ini,1)
                line_d = d_m_best*wall_d((k+1)/2);
                line_y_ = atan(j/line_d);
                line_y_ = (pi/2-line_y_)*im_h/pi;
                cor_opt = [cor_opt;cor_ini(k,1) line_y_];
            end
            % compute score
            score = interp2(corn,cor_opt(:,1),cor_opt(:,2));
            score = sum(log(score));
            if 1
            % add ceiling score
        	cor_all = [cor_opt(1,:);cor_opt(2,:);cor_opt(2,:);cor_opt(3,:);
               	cor_opt(3,:);cor_opt(4,:);cor_opt(4,:);cor_opt(1,:)];
    		[ uv ] = coords2uv( cor_all, im_w, im_h );
    		[ xyz ] = uv2xyzN( uv);
    		[ lines ] = lineFromTwoPoint( xyz(1:2:end,:), xyz(2:2:end,:) );
    		im = zeros(im_h, im_w, 3);
    		[ panoEdgeC ] = paintParameterLine_my(lines, im_w, im_h, im);
    		score_cl = zeros(1,4);

    		score_cl(1) = max(max(edg2(:,:,2).*(panoEdgeC(:,:,1) == 1)));
    		score_cl(2) = max(max(edg2(:,:,2).*(panoEdgeC(:,:,1) == 2)));
    		score_cl(3) = max(max(edg2(:,:,2).*(panoEdgeC(:,:,1) == 3)));
    		score_cl(4) = max(max(edg2(:,:,2).*(panoEdgeC(:,:,1) == 4)));
    		score_cl = sum(log(score_cl));
    		score = score + 0.5*score_cl;
    	end

            if score > score_best
                score_best = score;
                cor_best = cor_opt;
                h_m_best = j;
            end
        end
        if score_best > score_tp
            cor_ini(1:2:end,:) = cor_best;
        else
            score_best = score_tp;
            h_m_best = h_m;
        end
        score_all = score_all + score_best;
end
