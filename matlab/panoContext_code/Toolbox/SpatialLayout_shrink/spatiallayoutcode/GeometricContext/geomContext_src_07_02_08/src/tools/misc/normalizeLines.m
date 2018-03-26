function lines = normalizeLines(lines, imsize)
% Normalizes line parameters according to image height.
% lines([x1 x2 y1 y2 angle r])

if size(lines, 1) > 0
    lines(:, 1:2) = (lines(:, 1:2)-imsize(2)/2) / imsize(1);%meansize;
    lines(:, 3:4) = (lines(:, 3:4) - imsize(1)/2) / imsize(1);%meansize;
    %lines(:, 6) = lines(:, 6) / meansize;
    lines(:, 6) = mean(lines(:, 1:2), 2).*cos(lines(:, 5))+mean(lines(:, 3:4), 2).*sin(lines(:, 5));
end
