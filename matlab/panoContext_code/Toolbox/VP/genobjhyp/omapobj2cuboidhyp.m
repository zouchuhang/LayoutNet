function cuboidhyp = omapobj2cuboidhyp(rpo, reglabel, regori, vp, OMAP_FACTOR)

vps{1} = OMAP_FACTOR * vp{1};
vps{2} = OMAP_FACTOR * vp{2};
vps{3} = OMAP_FACTOR * vp{3};

regoriobj = regori;
statsobj = regionprops(reglabel,'PixelIdxList','PixelList');

%%
rm = regionmask(vps, size(reglabel));
sr = segregion(rm, statsobj);

convexpair = find_convex_pair(rpo, sr, regoriobj);

%% get all bounding quadrilaterals of convex pairs
[r c v] = find(convexpair);

rl = unique([r; c]);
for i = 1:length(rl)
    boundingquad(rl(i)) = fit_bounding_quad(statsobj(rl(i)).PixelList, regori(rl(i)), vps);
end

% % display quadrilaterals
% figure;
% img = zeros(size(reglabel));
% for i = 1:length(rl)
%     img(statsobj(rl(i)).PixelIdxList) = 1;
%     imshow(img,[]);
%     hold on;
%     for j = 1:4
%     plot(boundingquad(rl(i)).pt(j,1), boundingquad(rl(i)).pt(j,2), ...
%         'x', 'Markersize', 10, 'Color', [1 0 0]);
%     text(boundingquad(rl(i)).pt(j,1), boundingquad(rl(i)).pt(j,2), ...
%         num2str(j), 'FontSize', 20, 'Color', [1 0 0]);
%     end
%     pause;
%     hold off
%     img(statsobj(rl(i)).PixelIdxList) = 0;
% end

if ~exist('boundingquad', 'var')
    cuboidhyp = [];
    return;
end

%%

partcubhyp = convexpair2partialcuboidhyp(convexpair, boundingquad);

% upscale 
for i = 1:length(partcubhyp)
    for j = 1:8
        if ~isempty(partcubhyp(i).junc3(j).pt)
            partcubhyp(i).junc3(j).pt = partcubhyp(i).junc3(j).pt / OMAP_FACTOR;
        end
    end
end

%% get cuboidhyp from partial cuboids
for i = 1:length(partcubhyp)
    prm = cuboid_pts2prm_incomp(partcubhyp(i), vp);
    pts = cuboid_prm2pts(prm, vp);
    cuboidhyp(i) = cuboid_pts2hyp(pts, vp);
end


%% display...
% global omapmore;
% % figure;
% % img = zeros(size(reglabel));
% [row col] = find(convexpair);
% for i = 1:length(row)
%     
%     rimg = disp_region(omapmore, reglabel, [row(i) col(i)], OMAP_FACTOR);
% 
%     for j = 1:8
%         figure; imshow(rimg); hold on;
% 
%         pt = cat(1, boundingquad(row(i)).pt);
%         for k=1:4, text(pt(k,1),pt(k,2),num2str(k), 'Color',[0 1 1]); end
%         pt = pt([1 2 4 3 1], :) / OMAP_FACTOR;
%         plot(pt(:,1), pt(:,2), 'm', 'LineWidth',2);
%         
%         pt = cat(1, boundingquad(col(i)).pt);
%         for k=1:4, text(pt(k,1),pt(k,2),num2str(k), 'Color',[0 1 1]); end
%         pt = pt([1 2 4 3 1], :) / OMAP_FACTOR;
%         plot(pt(:,1), pt(:,2), 'm', 'LineWidth',2);    
% 
% 
%         disp_cubes(cuboidhyp((i-1)*8+j), []);
%         pause;
%         close;
%     end
% end

%% discard invalid cuboids
valid = check_validity(cuboidhyp, vp);
cuboidhyp = cuboidhyp(valid);

%% drop thin cuboids
cuboidhyp = drop_thin_cuboids(cuboidhyp, vp);

%% this does not work yet.
% ch = convexpair2cuboidhyp(convexpair, statsobj, vp);




%%
function convexpair = find_convex_pair(rel, sr, regoriobj)
% rel{1}(i,j) = 1 if i is up, j is down
% rel{2}(i,j) = 1 if i is left, j is right
% rel{3}(i,j) = 1 if i is front, j is back
% sr[n x 4]: sr(i,j)=1 if segment i is in region j
% regoriobj(i): orientation of region
% convexpair = 1,2,3,4, or 5 depending on configuration of convex pair
%    = 1 if top & front & region3or4
%    = 2 if top & left & region4
%    = 3 if top & right & region3
%    = 4 if left & front & region4
%    = 5 if right & front & region3

convexpair = zeros(size(rel{1}));

[r c] = find(rel{1} | rel{2} | rel{3});
row = [r; c];
col = [c; r];

for i = 1:length(row)
    si = row(i);
    sj = col(i);
    orii = regoriobj(si);
    orij = regoriobj(sj);
    regioni = sr(si,:);
    regionj = sr(sj,:);
    if orii==1 && orij==3 && ... % top & front & region3or4
       rel{1}(si,sj)==1 && rel{3}(sj,si)==1 && ...
       (regioni(3)==1 || regioni(4)==1) && ...
       (regionj(3)==1 || regionj(4)==1) && ...
       (regioni(1)==0 && regioni(2)==0)
           convexpair(si,sj) = 1;
    elseif orii==1 && orij==2 && ... % top & left & region4
        rel{1}(si,sj)==1 && rel{2}(sj,si)==1 && ...
        regioni(4)==1 && regionj(4)==1 && ...
        regioni(1)==0 && regioni(2)==0
            convexpair(si,sj) = 2;
    elseif orii==1 && orij==2 && ... % top & right & region3
        rel{1}(si,sj)==1 && rel{2}(si,sj)==1 && ...
        regioni(3)==1 && regionj(3)==1 && ...
        regioni(1)==0 && regioni(2)==0
            convexpair(si,sj) = 3;
    elseif orii==2 && orij==3 && ... % left & front & region4
        rel{2}(si,sj)==1 && rel{3}(sj,si)==1 && ...
        regioni(4)==1 && regionj(4)==1
            convexpair(si,sj) = 4;
    elseif orii==2 && orij==3 && ... & right & front & region3
        rel{2}(sj,si)==1 && rel{3}(sj,si)==1 && ...
        regioni(3)==1 && regionj(3)==1
            convexpair(si,sj) = 5;
    end
end

%%
function sr = segregion(rm, stats)
% sr[n x 4]: sr(i,j)=1 if segment i is in region j

sr = false(length(stats), 4);
for i = 1:length(stats)
    sr(i,1) = any( rm{1}(stats(i).PixelIdxList) );
    sr(i,2) = any( rm{2}(stats(i).PixelIdxList) );
    sr(i,3) = any( rm{3}(stats(i).PixelIdxList) );
    sr(i,4) = any( rm{4}(stats(i).PixelIdxList) );
end

%%
function rel = convert_rpo_to_udlrfb(rpo, vp)
% rpo{1,2,3}(j,k) = 1 if region j is closer to camera than k (k is closer to vp than j)
% rel{1}(i,j) = 1 if i is up, j is down
% rel{2}(i,j) = 1 if i is left, j is right
% rel{3}(i,j) = 1 if i is front, j is back

rel = rpo;
if vp{1}(2) < 0
    rel{1} = rel{1}';
end
if vp{2}(1) < 0
    rel{2} = rel{2}';
end


%%
function rm = regionmask(vp, imgsize)
% %    1 | 2
% %   ---+---
% %    3 | 4
% ***** TODO: fix this. need to check various conditions of vp *****

if vp{3}(2) < 1
    upmask = false(imgsize(1),imgsize(2));
elseif vp{3}(2) > imgsize(1)
    upmask = true(imgsize(1),imgsize(2));
else
    [p1ext p2ext] = extline(vp{2}, vp{3}, imgsize(2), imgsize(1));
    
    if p1ext(1) < p2ext(1)
        x = [p1ext(1) p2ext(1) imgsize(2)+100 -100];
    else
        x = [p1ext(1) p2ext(1) -100 imgsize(2)+100];
    end
    y = [p1ext(2) p2ext(2) -100 -100];
    
    upmask = poly2mask(x,y,imgsize(1),imgsize(2));
end
if vp{3}(1) < 1
    leftmask = false(imgsize(1),imgsize(2));
elseif vp{3}(1) > imgsize(2)
    leftmask = true(imgsize(1),imgsize(2));
else
    [p1ext p2ext] = extline(vp{1}, vp{3}, imgsize(2), imgsize(1));
    if p1ext(2) < p2ext(2)
        y = [p1ext(2) p2ext(2) imgsize(1)+100 -100];
    else
        y = [p1ext(2) p2ext(2) -100 imgsize(1)+100];
    end
    x = [p1ext(1) p2ext(1) -100 -100];
    
    leftmask = poly2mask(x,y,imgsize(1),imgsize(2));
end

rm{1} =  upmask &  leftmask;
rm{2} =  upmask & ~leftmask;
rm{3} = ~upmask &  leftmask;
rm{4} = ~upmask & ~leftmask;


