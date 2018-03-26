function cor_fn = samplingPanoBox(cor_id, corn, edg, edg2, im_h, im_w, options)

    corn = im2double(corn);
    edg = im2double(edg);

    score_fin = -Inf;
    score_fin_t = -Inf;
    search_grid = 5;
    line = cor_id(1:2:end,1)';
    line_can = repmat(line, size(line,2)*(2*search_grid+1),1);
    score_ln = interp2(corn,cor_id(:,1),cor_id(:,2));
    score_ln = log(score_ln);
    score_ln = [score_ln(1) + score_ln(2); score_ln(3) + score_ln(4);...
               score_ln(5) + score_ln(6); score_ln(7) + score_ln(8)];
    [~,score_id] = sort(score_ln);

    line_can(1:(2*search_grid+1),score_id(1)) = line(score_id(1))-search_grid:line(score_id(1))+search_grid;
    line_can((2*search_grid+1)+1:(2*search_grid+1)*2, score_id(2)) = line(score_id(2))-search_grid:line(score_id(2))+search_grid;
    line_can((2*search_grid+1)*2+1:(2*search_grid+1)*3, score_id(3)) = line(score_id(3))-search_grid:line(score_id(3))+search_grid;
    line_can((2*search_grid+1)*3+1:(2*search_grid+1)*4, score_id(4)) = line(score_id(4))-search_grid:line(score_id(4))+search_grid;
    edg2 = im2double(edg2);
    
    % no sampling
    %line_can = line;
    for line_n = 1:size(line_can,1)
        %disp(line_n)
        cor_ini = cor_id;
        cor_ini(1:2:end,1) = line_can(line_n,:);
        cor_ini(2:2:end,1) = line_can(line_n,:);
        [cor_ini, score_all] = sample_opt_pano_joint(cor_ini, im_w, im_h, corn, edg, edg2, options);
        if score_all > score_fin_t
            score_fin_t = score_all;
            cor_fn_t = cor_ini;
        end
        if mod(line_n,search_grid*2+1) == 0
            line_id = line_n/(search_grid*2+1);
            line_can(line_n+1:end,score_id(line_id)) = cor_fn_t(2*score_id(line_id),1);
            if score_fin_t > score_fin
                score_fin = score_fin_t;
                cor_fn = cor_fn_t;
            end
            score_fin_t = -Inf;
        end
    end
    % no sampling
    %if score_fin_t > score_fin
    %    score_fin = score_fin_t;
    %    cor_fn = cor_fn_t;
    %end
end