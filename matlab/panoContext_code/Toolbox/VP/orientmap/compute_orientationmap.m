function [omap omapstrict1 omapstrict2] = compute_orientationmap(lines, vp, imgsize, OMAP_FACTOR)


%% resize lines
lines_rsz = lines;
for i = 1:length(lines)
    lines_rsz(i).point1 = lines_rsz(i).point1 * OMAP_FACTOR;
    lines_rsz(i).point2 = lines_rsz(i).point2 * OMAP_FACTOR;
    lines_rsz(i).lineeq = lines_rsz(i).lineeq * OMAP_FACTOR;
end
vp{1} = vp{1} * OMAP_FACTOR;
vp{2} = vp{2} * OMAP_FACTOR;
vp{3} = vp{3} * OMAP_FACTOR;
omapsize = ceil(imgsize * OMAP_FACTOR);

%%
lineextimg = orient_from_lines(lines_rsz, vp, omapsize(2), omapsize(1));

%%
ao23 = (lineextimg{2,3,1} | lineextimg{2,3,2});
ao32 = (lineextimg{3,2,1} | lineextimg{3,2,2});
ao13 = (lineextimg{1,3,1} | lineextimg{1,3,2});
ao31 = (lineextimg{3,1,1} | lineextimg{3,1,2});
ao12 = (lineextimg{1,2,1} | lineextimg{1,2,2});
ao21 = (lineextimg{2,1,1} | lineextimg{2,1,2});
aa23 = (lineextimg{2,3,1} & lineextimg{2,3,2});
aa32 = (lineextimg{3,2,1} & lineextimg{3,2,2});
aa13 = (lineextimg{1,3,1} & lineextimg{1,3,2});
aa31 = (lineextimg{3,1,1} & lineextimg{3,1,2});
aa12 = (lineextimg{1,2,1} & lineextimg{1,2,2});
aa21 = (lineextimg{2,1,1} & lineextimg{2,1,2});

%% regular
% a{1} = (lineextimg{2,3,1} | lineextimg{2,3,2}) & (lineextimg{3,2,1} | lineextimg{3,2,2});
% a{2} = (lineextimg{1,3,1} | lineextimg{1,3,2}) & (lineextimg{3,1,1} | lineextimg{3,1,2});
% a{3} = (lineextimg{1,2,1} | lineextimg{1,2,2}) & (lineextimg{2,1,1} | lineextimg{2,1,2});
a{1} = ao23 & ao32;
a{2} = ao13 & ao31;
a{3} = ao12 & ao21;

b{1} = a{1} & ~a{2} & ~a{3};
b{2} = ~a{1} & a{2} & ~a{3};
b{3} = ~a{1} & ~a{2} & a{3};

omap = cat(3, b{1}, b{2}, b{3});


%% AND thing
% a{1} = lineextimg{2,3,1} & lineextimg{2,3,2} & lineextimg{3,2,1} & lineextimg{3,2,2};
% a{2} = lineextimg{1,3,1} & lineextimg{1,3,2} & lineextimg{3,1,1} & lineextimg{3,1,2};
% a{3} = lineextimg{1,2,1} & lineextimg{1,2,2} & lineextimg{2,1,1} & lineextimg{2,1,2};
a{1} = aa23 & aa32;
a{2} = aa13 & aa31;
a{3} = aa12 & aa21;

b{1} = a{1} & ~a{2} & ~a{3};
b{2} = ~a{1} & a{2} & ~a{3};
b{3} = ~a{1} & ~a{2} & a{3};

omapstrict2 = cat(3, b{1}, b{2}, b{3});

%% between
a{1} = (ao23 & aa32) | (aa23 & ao32);
a{2} = (ao13 & aa31) | (aa13 & ao31);
a{3} = (ao12 & aa21) | (aa12 & ao21);

b{1} = a{1} & ~a{2} & ~a{3};
b{2} = ~a{1} & a{2} & ~a{3};
b{3} = ~a{1} & ~a{2} & a{3};

omapstrict1 = cat(3, b{1}, b{2}, b{3});

