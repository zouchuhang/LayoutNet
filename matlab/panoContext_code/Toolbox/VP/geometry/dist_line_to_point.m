function d = dist_line_to_point(linep1, linep2, p)
% compute distance of point to line
% line is defined by 2 points

n = norm(linep1 - linep2);

if n<1e-10
    warning('dist_line_to_point.m: degenerate input');
end

d = abs((linep1(1)-linep2(1))*(linep1(2)-p(2)) ...
    -   (linep1(2)-linep2(2))*(linep1(1)-p(1))) ./ n;

