function im_seg = getSegMask_eval(cor_id, im_h, im_w)
% produce layout depth and segmentation result from predicted 2d 
% corner position in the panorama
    addpath(genpath('./affine_fit/'));
    addpath(genpath('./plane_line_intersect/'));

    trans_eval;
    % recenter
    xyz = [cor_a_h;cor_a; cor_b_h;cor_b; cor_c_h;cor_c; cor_d_h;cor_d];
    xyz(:,2) = xyz(:,2) - c_h;
    xyz = [xyz(:,1) xyz(:,3) xyz(:,2)];
    [uv] = coords2uv(cor_id, im_w, im_h);
    % top plane
    [n1,~,p1] = affine_fit(xyz(1:2:end,:));
    % bottom plane
    [n2,~,p2] = affine_fit(xyz(2:2:end,:));
    % 1st vertical plane
    [n3,~,p3] = affine_fit(xyz([1,2,3,4],:));
    % 2nd vertical plane
    [n4,~,p4] = affine_fit(xyz([3,4,5,6],:));
    % 3nd vertical plane
    [n5,~,p5] = affine_fit(xyz([5,6,7,8],:));
    % 4th vertical plane
    [n6,~,p6] = affine_fit(xyz([7,8,1,2],:));

    im_seg = zeros(im_h, im_w);
    [im_X,im_Y] = meshgrid(1:im_w, 1:im_h);
    im_cor = [im_X(:),im_Y(:)];
    [uv_im] = coords2uv(im_cor, im_w, im_h);
    [ xyz_im ] = uv2xyzN(uv_im);
    cen = [0 0 0];
    xyz_im = bsxfun(@plus, xyz_im, cen);
    for j = 1:size(xyz_im,1)
        [I1,check1]=plane_line_intersect(n1',p1,cen,xyz_im(j,:)*1000);
        [I2,check2]=plane_line_intersect(n2',p2,cen,xyz_im(j,:)*1000);
        [I3,check3]=plane_line_intersect(n3',p3,cen,xyz_im(j,:)*1000);
        [I4,check4]=plane_line_intersect(n4',p4,cen,xyz_im(j,:)*1000);
        [I5,check5]=plane_line_intersect(n5',p5,cen,xyz_im(j,:)*1000);
        [I6,check6]=plane_line_intersect(n6',p6,cen,xyz_im(j,:)*1000);
        check = [check1, check2, check3, check4, check5, check6];
        I = [I1;I2;I3;I4;I5;I6];
        id = find(check == 1);
        dist = I(id,:);
        dist = sqrt(sum(dist.*dist,2));
        [dep, idx] = min(dist);
        im_seg(im_cor(j,2), im_cor(j,1)) = id(idx);
    end
    im_seg(im_seg>2) = 3;

end
