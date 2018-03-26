function im = draw_circle_image(im, cx, cy, r, color, varargin)

b =[0];
if length(varargin)>0
    w = varargin{1};
    b = [ceil(-w/2):ceil(w/2)-1];   
end

if numel(color)==3
    color = reshape(color, [1 1 3]);
end

theta = [0:1/(2*r):2*pi];
x = round(cos(theta)*r + cx);
y = round(sin(theta)*r + cy);

ind = find( (x>0) & (y>0) & (x<=size(im, 2)) & (y<=size(im, 1)));
x = x(ind);
y = y(ind);

for i = 1:length(x)    
    im(y(i)+b, x(i)+b, :) = repmat(color, [numel(b) numel(b) 1]);
end

