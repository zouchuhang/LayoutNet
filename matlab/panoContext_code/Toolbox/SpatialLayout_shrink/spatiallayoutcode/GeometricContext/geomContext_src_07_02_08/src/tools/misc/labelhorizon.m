function horizons = labelhorizon(imdir, fn)
f = length(horizons);
for f = 1:length(fn)
    im = imread([imdir '/' fn{f}]);
    
    ok = 0;
    while(ok~=121)
        hold off
        figure(1), imshow(im)
        [x y]= ginput(1);
        %m = (y(2) - y(1)) / (x(2) - x(1));
        m = 0;
        b = y(1) - m*x(1);
        disp('is this ok? (y/n)')
        hold on, plot([1 size(im, 2)], [m+b m*size(im, 2)+b]); 
        [x y ok] = ginput(1);
    end
    horizons(f).m = m;
    horizons(f).b = b;
    horizons(f).imname = fn{f};
end
    
    