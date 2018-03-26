function d = sqdist(x, y)
% L2 dist from x to y
% x: [1 x d]
% y: [n x d]
% d: [n x 1]

d = sqrt(sum((repmat(x, size(y,1), 1) - y).^2, 2));

