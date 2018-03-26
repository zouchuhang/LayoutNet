function [cubcompat closercub] = cuboid_cuboid_compat(cuboidhyp, imgsize, OMAP_FACTOR, vp)
% buildinghyp: length==1 must be 1
% cuboidhyp: length>=1 can be many
% cubcompat: cubcompat(i,j)==1 if cube i and j do not occupy same 3D space
% closercub: closercub(i,j)==1 if cube i is closer in 3D than cube j

num_cuboid = length(cuboidhyp);

omapsize = ceil(imgsize * OMAP_FACTOR);

scandir = 3;
scanlineendpts = vpscanlines(vp, scandir, omapsize, OMAP_FACTOR);


cubcompat = false(num_cuboid, num_cuboid);
closercub = sparse(num_cuboid, num_cuboid);
% cube 1 bottom
for i = 1:num_cuboid
    pt = cat(1, cuboidhyp(i).junc3([3 4 8 7]).pt);
    pt = pt * OMAP_FACTOR;
    cuboid_bottom_1 = poly2mask(pt(:,1), pt(:,2), omapsize(1), omapsize(2));

    % cube 2 bottom
    for j = i+1:num_cuboid
        pt = cat(1, cuboidhyp(j).junc3([3 4 8 7]).pt);
        pt = pt * OMAP_FACTOR;
        cuboid_bottom_2 = poly2mask(pt(:,1), pt(:,2), omapsize(1), omapsize(2));

        if any( cuboid_bottom_1(:)==1 & cuboid_bottom_2(:)==1 )
            % cubcompat(i,j) = false;
            % cubcompat(j,i) = false;
        else
            cubcompat(i,j) = true;
            cubcompat(j,i) = true;
            cc = closer_cuboid(cuboid_bottom_1, cuboid_bottom_2, scanlineendpts);
            if cc==1
                closercub(i,j) = 1; % i is closer than j
            elseif cc==2
                closercub(j,i) = 1; % j is closer than i
            end
        end
    end
end

%%
function closercub = closer_cuboid(cuboid_bottom_1, cuboid_bottom_2, scanlineendpts)
% closercub: 0 if can't tell
%            1 if cub1 is closer
%            2 if cub2 is closer

for i = 1:size(scanlineendpts,1)
    scanlinesamples = endpts2samples(scanlineendpts(i,:));
    
    cubscan1 = cuboid_bottom_1(sub2ind(size(cuboid_bottom_1),scanlinesamples(:,2),scanlinesamples(:,1)));
    cubscan2 = cuboid_bottom_2(sub2ind(size(cuboid_bottom_2),scanlinesamples(:,2),scanlinesamples(:,1)));

    min1 = find(cubscan1==1, 1, 'first');
    max1 = find(cubscan1==1, 1, 'last');
    min2 = find(cubscan2==1, 1, 'first');
    max2 = find(cubscan2==1, 1, 'last');
    
    if     max1<min2
        closercub = 1; % 1 is closer to camera than 2
        break;
    elseif max2<min1
        closercub = 2; % 2 is closer to camera than 1
        break;
    else
        closercub = 0;
    end
end


%%
function samples = endpts2samples(endpts)
d = norm(endpts(1:2) - endpts(3:4));
samples = [linspace(endpts(1),endpts(3),d)' ...
    linspace(endpts(2),endpts(4),d)'];
samples = round(samples);


