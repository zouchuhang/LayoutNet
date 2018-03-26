function image = draw_line_image2(image, lines, intensity, varargin)
% image(height, width, length(intensity))
% lines([x1 x2 y1 y2], num_lines)
% intensity([r g b]) or intensity([gray])

b = [0];
if length(varargin)>0
    w = varargin{1};  
    b = [ceil(-w/2):ceil(w/2)-1];  
%    b =[ceil(-w/2):floor(w/2)+w-1]; 
%    disp(num2str(b))
end
maxb = 0;
minb = 0;

intensity = reshape(intensity, [1 1 numel(intensity)]);

for index = 1:size(lines,2)
    %x1 = min(max(round(lines(1, index)), 1-minb), size(image, 2)-maxb);
    %x2 = min(max(round(lines(2, index)), 1-minb), size(image, 2)-max);
    %y1 = min(max(round(lines(3, index)), 1-minb), size(image, 1)-maxb);
    %y2 = min(max(round(lines(4, index)), 1-minb), size(image, 1)-maxb);       
    x1 = round(lines(1, index));
    x2 = round(lines(2, index));
    y1 = round(lines(3, index));
    y2 = round(lines(4, index));
    
    line_width = abs(x2-x1);
    line_height = abs(y2-y1);

    if line_width >= line_height
        for x = x1:sign(x2-x1):x2
            y = round((x-x1)/(x2-x1)*(y2-y1)+y1);
            try
                image(y+b, x, :) = repmat(intensity, [numel(b) 1 1]);
            catch
            end
        end
    else
        for y = y1:sign(y2-y1):y2
            x = round((y-y1)/(y2-y1)*(x2-x1)+x1);
            try
                image(y, x+b, :) = repmat(intensity, [1 numel(b) 1]);
            catch
            end
        end
    end
end