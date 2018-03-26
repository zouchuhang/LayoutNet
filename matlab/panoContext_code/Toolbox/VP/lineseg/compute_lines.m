function [lines linesmore] = compute_lines(img)

diag = sqrt(size(img,1).^2 + size(img,2).^2);

longlen = ceil(diag/30);
[edgeimg lines] = canny_wrapper(img, 5, 10, longlen, 1);

shortlen = ceil(diag/50);
[edgeimg linesmore] = canny_wrapper(img, 2, 5, shortlen, 1);

for i=1:length(lines)
    lines(i).length = norm(lines(i).point1-lines(i).point2);
end
