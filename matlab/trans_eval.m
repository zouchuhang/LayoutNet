c_h = 1;
%im_h = 512;
%im_w = 1024;

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

    if 0
    cor_e_x = cor_id(10,1) - im_w/2;
    cor_e_y = (im_h - cor_id(10,2))-im_h/2;
    cor_e_y_ = (im_h - cor_id(9,2))-im_h/2;
    theta_x = 2*pi*cor_e_x/im_w;
    theta_y = pi*cor_e_y/im_h;
    theta_y_ = pi*cor_e_y_/im_h;
    r = abs(cot(theta_y))*c_h;
    cor_e_X = r*cos(-theta_x+pi/2);% *ch
    cor_e_Z = r*sin(-theta_x+pi/2);% *ch
    cor_e_Y_ = r*tan(theta_y_) +c_h;

    cor_f_x = cor_id(12,1) - im_w/2;
    cor_f_y = (im_h - cor_id(12,2))-im_h/2;
    cor_f_y_ = (im_h - cor_id(11,2))-im_h/2;
    theta_x = 2*pi*cor_f_x/im_w;
    theta_y = pi*cor_f_y/im_h;
    theta_y_ = pi*cor_f_y_/im_h;
    r = abs(cot(theta_y))*c_h;
    cor_f_X = r*cos(-theta_x+pi/2);% *ch
    cor_f_Z = r*sin(-theta_x+pi/2);% *ch
    cor_f_Y_ = r*tan(theta_y_) +c_h;
    end
    
    cor_a = [cor_a_X, 0, cor_a_Z];
    cor_b = [cor_b_X, 0, cor_b_Z];
    cor_c = [cor_c_X, 0, cor_c_Z];
    cor_d = [cor_d_X, 0, cor_d_Z];
    %cor_e = [cor_e_X, 0, cor_e_Z];
    %cor_f = [cor_f_X, 0, cor_f_Z];

    box_h = mean([cor_a_Y_, cor_b_Y_, cor_c_Y_, cor_d_Y_]);%, cor_e_Y_, cor_f_Y_]);

    cor_a_h = [cor_a_X, box_h, cor_a_Z];
    cor_b_h = [cor_b_X, box_h, cor_b_Z];
    cor_c_h = [cor_c_X, box_h, cor_c_Z];
    cor_d_h = [cor_d_X, box_h, cor_d_Z];
    %cor_e_h = [cor_e_X, box_h, cor_e_Z];
    %cor_f_h = [cor_f_X, box_h, cor_f_Z];
