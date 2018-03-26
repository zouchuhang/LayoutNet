function image = draw_arrow_image(image, arrows, intensity, varargin)
% image(height, width)
% arrow.{x, y, angle (in deg), radius, head_length, head_base_angle}
% intensity([r g b]) or intensity([gray])

if isempty(varargin)
    varargin{1} = 1;
end

for a = 1:length(arrows)

    arrows(a).angle = arrows(a).angle / 180 * pi;
    arrows(a).head_base_angle = arrows(a).head_base_angle / 180 * pi;

    % form shaft of arrow
    x1 = arrows(a).x;
    y1 = arrows(a).y;    
    x2 = x1+cos(arrows(a).angle)*arrows(a).radius;
    y2 = y1-sin(arrows(a).angle)*arrows(a).radius;    
    lines(1).x = [x1 x2];
    lines(1).y = [y1 y2];    
        
    % form head of arrow
    x1 = x2;  % start at tip
    y1 = y2;
    % left side of head    
    angle = arrows(a).angle - pi - arrows(a).head_base_angle;
    x2 = x1+cos(angle)*arrows(a).head_length;
    y2 = y1-sin(angle)*arrows(a).head_length;
    lines(2).x = [x1 x2];
    lines(2).y = [y1 y2]; 
    % right side of head  
    angle = arrows(a).angle - pi + arrows(a).head_base_angle;
    x2 = x1+cos(angle)*arrows(a).head_length;
    y2 = y1-sin(angle)*arrows(a).head_length;
    lines(3).x = [x1 x2];
    lines(3).y = [y1 y2]; 
    
    % draw lines
    image = draw_line_image(image, lines, intensity, varargin{1});
    
end