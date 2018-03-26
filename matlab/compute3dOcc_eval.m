function score = compute3dOcc_eval(gt, cor_id, im_h, im_w);
% compute 3d IoU between gt corner position and 
% predicted corner posision in 2D

    addpath(genpath('./mesh2voxel/'));

    c_h = 1.7;
    voxelScale = 200;

    % gen box gt
    cor_a_x = cor_id(2,1) - im_w/2;
    cor_a_y = (im_h - cor_id(2,2))-im_h/2;
    cor_a_y_ = (im_h - cor_id(1,2))-im_h/2;
    theta_x = 2*pi*cor_a_x/im_w;
    theta_y = pi*cor_a_y/im_h;
    theta_y_ = pi*cor_a_y_/im_h;
    r = abs(cot(theta_y))*c_h;
    cor_a_X = r*cos(-theta_x+pi/2);% *ch
    cor_a_Z = r*sin(-theta_x+pi/2);% *ch 
    cor_a_Y_ = r*tan(theta_y_)+c_h;

    cor_b_x = cor_id(4,1) - im_w/2;
    cor_b_y = (im_h - cor_id(4,2))-im_h/2;
    cor_b_y_ = (im_h - cor_id(3,2))-im_h/2;
    theta_x = 2*pi*cor_b_x/im_w;
    theta_y = pi*cor_b_y/im_h;
    theta_y_ = pi*cor_b_y_/im_h;
    r = abs(cot(theta_y))*c_h;
    cor_b_X = r*cos(-theta_x+pi/2);% *ch
    cor_b_Z = r*sin(-theta_x+pi/2);% *ch 
    cor_b_Y_ = r*tan(theta_y_) +c_h;

    cor_c_x = cor_id(6,1) - im_w/2;
    cor_c_y = (im_h - cor_id(6,2))-im_h/2;
    cor_c_y_ = (im_h - cor_id(5,2))-im_h/2;
    theta_x = 2*pi*cor_c_x/im_w;
    theta_y = pi*cor_c_y/im_h;
    theta_y_ = pi*cor_c_y_/im_h;
    r = abs(cot(theta_y))*c_h;
    cor_c_X = r*cos(-theta_x+pi/2);% *ch
    cor_c_Z = r*sin(-theta_x+pi/2);% *ch 
    cor_c_Y_ = r*tan(theta_y_) +c_h;

    cor_d_x = cor_id(8,1) - im_w/2;
    cor_d_y = (im_h - cor_id(8,2))-im_h/2;
    cor_d_y_ = (im_h - cor_id(7,2))-im_h/2;
    theta_x = 2*pi*cor_d_x/im_w;
    theta_y = pi*cor_d_y/im_h;
    theta_y_ = pi*cor_d_y_/im_h;
    r = abs(cot(theta_y))*c_h;
    cor_d_X = r*cos(-theta_x+pi/2);% *ch
    cor_d_Z = r*sin(-theta_x+pi/2);% *ch
    cor_d_Y_ = r*tan(theta_y_) +c_h;

    % check correctness
    cor_a = [cor_a_X, 0, cor_a_Z];
    cor_b = [cor_b_X, 0, cor_b_Z];
    cor_c = [cor_c_X, 0, cor_c_Z];
    cor_d = [cor_d_X, 0, cor_d_Z];
   
    cor_a_h = [cor_a_X, cor_a_Y_, cor_a_Z];
    cor_b_h = [cor_b_X, cor_b_Y_, cor_b_Z];
    cor_c_h = [cor_c_X, cor_c_Y_, cor_c_Z];
    cor_d_h = [cor_d_X, cor_d_Y_, cor_d_Z];

    
    % gt part
    % gen box gt
    gt_a_x = gt(2,1) - im_w/2;
    gt_a_y = (im_h - gt(2,2))-im_h/2;
    gt_a_y_ = (im_h - gt(1,2))-im_h/2;
    theta_x = 2*pi*gt_a_x/im_w;
    theta_y = pi*gt_a_y/im_h;
    theta_y_ = pi*gt_a_y_/im_h;
    r = abs(cot(theta_y))*c_h;
    gt_a_X = r*cos(-theta_x+pi/2);% *ch
    gt_a_Z = r*sin(-theta_x+pi/2);% *ch 
    gt_a_Y_ = r*tan(theta_y_)+c_h;

    gt_b_x = gt(4,1) - im_w/2;
    gt_b_y = (im_h - gt(4,2))-im_h/2;
    gt_b_y_ = (im_h - gt(3,2))-im_h/2;
    theta_x = 2*pi*gt_b_x/im_w;
    theta_y = pi*gt_b_y/im_h;
    theta_y_ = pi*gt_b_y_/im_h;
    r = abs(cot(theta_y))*c_h;
    gt_b_X = r*cos(-theta_x+pi/2);% *ch
    gt_b_Z = r*sin(-theta_x+pi/2);% *ch 
    gt_b_Y_ = r*tan(theta_y_) +c_h;

    gt_c_x = gt(6,1) - im_w/2;
    gt_c_y = (im_h - gt(6,2))-im_h/2;
    gt_c_y_ = (im_h - gt(5,2))-im_h/2;
    theta_x = 2*pi*gt_c_x/im_w;
    theta_y = pi*gt_c_y/im_h;
    theta_y_ = pi*gt_c_y_/im_h;
    r = abs(cot(theta_y))*c_h;
    gt_c_X = r*cos(-theta_x+pi/2);% *ch
    gt_c_Z = r*sin(-theta_x+pi/2);% *ch 
    gt_c_Y_ = r*tan(theta_y_) +c_h;

    gt_d_x = gt(8,1) - im_w/2;
    gt_d_y = (im_h - gt(8,2))-im_h/2;
    gt_d_y_ = (im_h - gt(7,2))-im_h/2;
    theta_x = 2*pi*gt_d_x/im_w;
    theta_y = pi*gt_d_y/im_h;
    theta_y_ = pi*gt_d_y_/im_h;
    r = abs(cot(theta_y))*c_h;
    gt_d_X = r*cos(-theta_x+pi/2);% *ch
    gt_d_Z = r*sin(-theta_x+pi/2);% *ch
    gt_d_Y_ = r*tan(theta_y_) +c_h;

    gt_a = [gt_a_X, 0, gt_a_Z];
    gt_b = [gt_b_X, 0, gt_b_Z];
    gt_c = [gt_c_X, 0, gt_c_Z];
    gt_d = [gt_d_X, 0, gt_d_Z];
    
    gt_a_h = [gt_a_X, gt_a_Y_, gt_a_Z];
    gt_b_h = [gt_b_X, gt_b_Y_, gt_b_Z];
    gt_c_h = [gt_c_X, gt_c_Y_, gt_c_Z];
    gt_d_h = [gt_d_X, gt_d_Y_, gt_d_Z];


    % box
    fvc.vertices = [cor_a; cor_b; cor_c; cor_d;cor_a_h; cor_b_h;cor_c_h;cor_d_h];
    fvc.faces = [1 2 3; 1 3 4;
                 5 6 7; 5 7 8;
                 1 2 5; 2 5 7;
                 3 4 8; 3 8 7;
                 2 3 7; 2 7 6;
                 1 4 5; 4 5 8];


    % gt
    fvc_gt.vertices = [gt_a; gt_b; gt_c; gt_d;gt_a_h; gt_b_h;gt_c_h;gt_d_h];
    fvc_gt.faces = [1 2 3; 1 3 4;
                 5 6 7; 5 7 8;
                 1 2 5; 2 5 7;
                 3 4 8; 3 8 7;
                 2 3 7; 2 7 6;
                 1 4 5; 4 5 8];


    fmax = max([fvc.vertices; fvc_gt.vertices]);
    fmin = min([fvc.vertices; fvc_gt.vertices]);
    
    dist = pdist2([fvc.vertices; fvc_gt.vertices], [fvc.vertices; fvc_gt.vertices]);
    [dist,~] = max(dist(:));
    dist = sqrt(dist) * 1.1;

    fvc.vertices = bsxfun(@minus, fvc.vertices,(fmax+fmin)/2);
    fvc.vertices = fvc.vertices*voxelScale/dist;
    pred_vox = polygon2voxel(fvc, [1 1 1]*voxelScale, 'center');


    fvc_gt.vertices = bsxfun(@minus, fvc_gt.vertices,(fmax+fmin)/2);
    fvc_gt.vertices = fvc_gt.vertices*voxelScale/dist;
    gt_vox = polygon2voxel(fvc_gt, [1 1 1]*voxelScale, 'center');

    [ score ] = VoxelIOUScore( pred_vox, gt_vox );

end
