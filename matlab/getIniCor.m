function [cor_id, pks, pk_loc] = getIniCor(cor_m, corn, edg_m, im_h)

    [pks_c,locs_c,w_c,p_c] = findpeaks(cor_m,'MinPeakProminence',58,'MinPeakDistance',20);
    [~, pk_id_c] = sort(pks_c, 'descend');
    pk_loc_c = locs_c(pk_id_c(1:min(4,numel(pks_c))));
    pk_loc_c = sort(pk_loc_c);
    if numel(pk_loc_c) < 4
        pk_loc_c = [1 pk_loc_c];
    end
    [pks,locs,w,p] = findpeaks(edg_m,'MinPeakProminence',20,'MinPeakDistance',20);
    
    [~, pk_id] = sort(pks, 'descend');
    pk_loc = locs(pk_id(1:min(4,numel(pks))));
    pk_loc = sort(pk_loc);
    
    cor_id = [];
    for j = 1:4
       [pks_t,locs_t,w_t,p_t] = findpeaks(double(corn(:,pk_loc_c(j))),'MinPeakProminence',50,'MinPeakDistance',20);
       if numel(pks_t) < 2
           [pks_t,locs_t,w_t,p_t] = findpeaks(double(corn(:,pk_loc_c(j))),'MinPeakProminence',5,'MinPeakDistance',20);
       end
       if numel(pks_t) < 2
            locs_t = [im_h/2,im_h/2];
            pks_t = [0,0];
       end
       [~, pk_id_t] = sort(pks_t, 'descend');
       pk_loc_t = locs_t(pk_id_t(1:min(2,numel(pks_t))));
       pk_loc_t = sort(pk_loc_t);
       cor_id = [cor_id ; pk_loc_c(j) pk_loc_t(1); pk_loc_c(j) pk_loc_t(2)];
    end
end