function im = draw_x_image(im, xc, yc, s, color, varargin)

x1 = xc - s/2;
y1 = yc - s/2;
im = draw_line_image2(im, [x1 x1+s y1 y1+s; x1 x1+s y1+s y1]', color, varargin{1});
