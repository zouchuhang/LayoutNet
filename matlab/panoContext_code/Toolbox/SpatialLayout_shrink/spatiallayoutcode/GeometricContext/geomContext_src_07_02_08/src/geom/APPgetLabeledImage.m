function image = APPgetLabeledImage(image, imsegs, vLabels, vConf, hLabels, hConf)
% image = APPgetLabeledImage(image, imsegs, vLabels, vConf, hLabels, hConf)
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005
%image = image / double(max(image(:)));

scale = size(image, 1) / size(vLabels, 1);

image = rgb2hsv(image);
image(:, :, 2) = 0.0*image(:, :, 2);
image = hsv2rgb(image);
           
drawn = zeros(size(vLabels));

rad = round(max(size(image))/30);
markW = ceil(rad/10);
right_arrow = struct('x', 0, 'y', 0, 'angle', 0, 'radius', rad, 'head_length', rad/2, 'head_base_angle', 30);
up_arrow = struct('x', 0, 'y', 0, 'angle', 90, 'radius', rad, 'head_length', rad/2, 'head_base_angle', 30);
up_left_arrow = struct('x', 0, 'y', 0, 'angle', 135, 'radius', rad, 'head_length', rad/2, 'head_base_angle', 30);    
up_right_arrow = struct('x', 0, 'y', 0, 'angle', 45, 'radius', rad, 'head_length', rad/2, 'head_base_angle', 30);    
left_arrow = struct('x', 0, 'y', 0, 'angle', 180, 'radius', rad, 'head_length', rad/2, 'head_base_angle', 30);
down_arrow = struct('x', 0, 'y', 0, 'angle', 270, 'radius', rad, 'head_length', rad/2, 'head_base_angle', 30);     
toward_arrow = struct('x', 0, 'y', 0, 'angle', 250, 'radius', rad*2/3, 'head_length', rad/3, 'head_base_angle', 60);

nsegs = imsegs.nseg;
sinds = APPgetSpIndsOld((1:nsegs), imsegs);
for s = 1:nsegs

    ind = sinds(s).inds;
    
    sat = max(min(vConf(s), 1), 0.33);

    val = 1;
    mix_rat = 0.5;
    arrows = [];
    hue = 0;
    if strcmp(vLabels{s}, 'sky')
        hue = 177/255;
    elseif strcmp(vLabels{s}, '000')
        hue = 112/255;
        % elseif strcmp(vLabels{s}, 'rnd') | strcmp(hLabels{s}, 'rnd')
       % hue = 226/255;        
    elseif strcmp(vLabels{s}, '045')
        %hue = 45/255;
        hue = 0/255;
    elseif strcmp(vLabels{s}, '090')
        hue = 0/255;
        %hue = 0/255;
    elseif strcmp(vLabels{s}, 'dst')
        %hue = 80/255;
        %sat = sat/2;
        hue = 0/255;
    else
        disp(vLabels{s})
        sat = 0;
    end
    %if strcmp(hLabels{s}, 'rnd') & (strcmp(vLabels{s}, '090') | strcmp(vLabels{s}, '045'))
    %    hue = 0/255;
    %    sat = sat/2;
    %end
   
    
        
    width = size(image, 2);
    height = size(image, 1);

    % fill in region with color determined by label
    colors = hsv2rgb([hue sat val]);
    label_colors = colors;
    
    colors = reshape(colors, [1 1 3]);

    y_ind = mod(ind-1, height)+1;
    x_ind = floor((ind-1)/height)+1;
    for s = 1:length(ind)
        image(y_ind(s), x_ind(s), (1:3)) = (1-mix_rat)*image(y_ind(s), x_ind(s), :) + mix_rat*colors;
    end

end
if 1
    
for x = rad:round(rad*1.51):(width-rad)
    for y = rad:round(rad*1.51):(height-rad)
        
        hlabel = hLabels{imsegs.segimage(y, x)};
        vlabel = vLabels{imsegs.segimage(y, x)};
        
        intensity = 0;
        label_colors = [1 1 1];
        
        arrows = [];
        if strcmp(vlabel, '090') | strcmp(vlabel, '045')
            if strcmp(hlabel, '045')
                arrows = left_arrow;
            elseif strcmp(hlabel, '090')
                arrows = up_arrow;
            elseif strcmp(hlabel, '135')
                arrows = right_arrow;                
            elseif strcmp(hlabel, 'por')
                image = draw_circle_image(image, x, y, rad/2, 0, markW);
            elseif strcmp(hlabel, 'sol')
                image = draw_x_image(image, x, y, rad, 0, markW);
            end
                
        end
        
        % draw arrows                     
        for a = 1:length(arrows)
            arrows(a).x = x;
            arrows(a).y = y; 
        end

        image = draw_arrow_image(image, arrows, 0, markW);              
        
    end
     
end
end

