function image = draw_line_image(image, lines, intensity, varargin)
% image(height, width, length(intensity))
% lines(num_lines).{x([1 2]), y([1 2])}
% intensity([r g b]) or intensity([gray])

b = 0;
if length(varargin)>0
    w = varargin{1};
    b = [ceil(-w/2):ceil(w/2)-1];  
end

intensity = reshape(intensity, [1 1 numel(intensity)]);

[imh, imw, nb] = size(image);

for index = 1:length(lines)

%     x1 = min(max(round(lines(index).x(1)), 1), size(image, 2));
%     x2 = min(max(round(lines(index).x(2)), 1), size(image, 2));
%     y1 = min(max(round(lines(index).y(1)), 1), size(image, 1));
%     y2 = min(max(round(lines(index).y(2)), 1), size(image, 1));          
    
    x1 = round(lines(index).x(1));
    x2 = round(lines(index).x(2));
    y1 = round(lines(index).y(1));
    y2 = round(lines(index).y(2));     
    
    line_width = abs(x2-x1);
    line_height = abs(y2-y1);
    if line_width >= line_height        
        for x = min(max(x1,1),imw):sign(x2-x1):max(min(x2,imw),1)
            y = round((x-x1)/(x2-x1)*(y2-y1)+y1) + b;
            y = y((y>0) & (y<=imh));
            image(y, x, :) = repmat(intensity, [numel(y) 1 1]);
        end
    else
        for y = min(max(y1,1),imh):sign(y2-y1):max(min(y2, imh), 1)
            x = round((y-y1)/(y2-y1)*(x2-x1)+x1) + b;
            x = x((x>0) & (x<=imw));
            image(y, x, :) = repmat(intensity, [1 numel(x) 1]);
        end
    end
end

[imh2, imw2, nb] = size(image);
if ~all([imh imw]==[imh2 imw2])
    keyboard;
end
