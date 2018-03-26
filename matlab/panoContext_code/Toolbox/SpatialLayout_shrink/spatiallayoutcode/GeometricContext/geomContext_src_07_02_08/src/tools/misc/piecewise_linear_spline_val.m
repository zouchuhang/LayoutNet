function y = piecewise_linear_spline_val(x, params, s)
% gets the value from a piecewise linear spline
% params(nendpts, [m b x1 x2 maxy disconnected])
y = zeros(size(x));
ind = zeros(size(x));

if exist('s', 'var')
    ind = s;
else
    for i = 1:numel(x)      
        
        
        if x(i) <= params(1, 3)
            ind(i) = 1;
        elseif x(i) >= params(end, 4)
            ind(i) = size(params, 1);
        else
            tmp = find((x(i) >= params(:, 3)) & (x(i) <= params(:, 4)));
            ind(i) = tmp(1);                            
        end
    end
end
m = params(ind, 1);
b = params(ind, 2);

y = m.*x + b;