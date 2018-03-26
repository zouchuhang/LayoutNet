function image = draw_perimeter_image(image, map, intensity)
% draws a tight perimeter along ones in the map
% image(height, width, length(intensity))
% map(height, width, length(intensity))

scale = size(image, 1) / size(map, 1);
if size(image, 2) / size(map, 2) ~= scale
    disp('warning: map proportions not to scale!');
end

% get bounds
im_y = mod((0:size(map,2)*size(map, 1)-1), size(map, 1))+1;
im_x = floor((0:size(map,2)*size(map, 1)-1)/size(map, 1))+1;
solids_ind = find(map==1); 
solids_x = im_x(solids_ind);
solids_y = im_y(solids_ind);

min_x = min(solids_x);
min_y = min(solids_y);
max_x = max(solids_x);
max_y = max(solids_y);

% draw and bottom top of perimeter
for x = min_x:max_x
    if min_y == 1 & map(min_y, x)
        image = draw_line_image2(image, ([x (x+1) min_y min_y]-1)'*scale+1, intensity);
    end   
    if map(max_y, x)
            image = draw_line_image2(image, ([x (x+1) max_y+1 max_y+1]-1)'*scale+1, intensity);
    end       
    for y = max(min_y,2):max_y
      
        if map(y, x)==1 & map(y-1, x)~=1
            image = draw_line_image2(image, ([x (x+1) y y]-1)'*scale+1, intensity);
        end
        if map(y, x)~=1 & map(y-1, x)==1
            image = draw_line_image2(image, ([x (x+1) y y]-1)'*scale+1, intensity);
        end
    end
end

% draw left and right of perimeter
for y = min_y:max_y
    if min_x == 1 & map(y, min_x)
        image = draw_line_image2(image, ([min_x min_x y (y+1)]-1)'*scale+1, intensity);
    end    
    if map(y, max_x)
        image = draw_line_image2(image, ([max_x+1 max_x+1 y (y+1)]-1)'*scale+1, intensity);
    end     
    for x = max(min_x, 2):max_x        
        if map(y, x)==1 & map(y, x-1)~=1
            image = draw_line_image2(image, ([x x y (y+1)]-1)'*scale+1, intensity);
        end
        if map(y, x)~=1 & map(y, x-1)==1
            image = draw_line_image2(image, ([x x y (y+1)]-1)'*scale+1, intensity);
        end
    end
end



            
