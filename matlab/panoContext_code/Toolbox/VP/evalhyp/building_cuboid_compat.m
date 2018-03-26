function compat = building_cuboid_compat(buildinghyp, cuboidhyp, imgsize, vp, OMAP_FACTOR)
% buildinghyp: length==1 must be 1
% cuboidhyp: length>=1 can be many

num_cuboid = length(cuboidhyp);
num_building = length(buildinghyp);

compat = false(num_building, num_cuboid);

imgpoly = [1 1; imgsize(2) 1; imgsize(2) imgsize(1); 1 imgsize(1); 1 1];
[imgpoly(:,1) imgpoly(:,2)]= poly2cw(imgpoly(:,1), imgpoly(:,2));

omapsize = ceil(imgsize * OMAP_FACTOR);


for i = 1:num_building
%     [orientlabel_img regionlabel_img] = ...
%         get_labelimg_frombox(buildinghyp(i).box, imgsize(2), imgsize(1), vp);
%     building_bottom = (regionlabel_img==length(buildinghyp(i).box)+1);

    bp = building_floor_poly(buildinghyp(i), imgsize);
    
    bp = bp * OMAP_FACTOR;
    bld_bottom = poly2mask(bp(:,1), bp(:,2), omapsize(1), omapsize(2));
    
    % cuboid bottom
    for j = 1:num_cuboid
        cp = cat(1, cuboidhyp(j).junc3([4 3 7 8]).pt);
        
        cp = cp * OMAP_FACTOR;
        cub_bottom = poly2mask(cp(:,1), cp(:,2), omapsize(1), omapsize(2));
        cubearea = length(find(cub_bottom));
        cubebldintersect = length(find(bld_bottom==1 & cub_bottom==1));
        if cubebldintersect / cubearea > 0.95
            compat(i,j) = true;
        end
    end
end


% % % % % % % INTERSECTING POLYGONS AREN'T FASTER THAN MASKING...
% % % % % % % building bottom
% % % % % % for i = 1:num_building
% % % % % % %     [orientlabel_img regionlabel_img] = ...
% % % % % % %         get_labelimg_frombox(buildinghyp(i).box, imgsize(2), imgsize(1), vp);
% % % % % % %     building_bottom = (regionlabel_img==length(buildinghyp(i).box)+1);
% % % % % % 
% % % % % %     bp = building_floor_poly(buildinghyp(i), imgsize);
% % % % % % %     assert(ispolycw(bp(:,1),bp(:,2)));
% % % % % %     [bp(:,1) bp(:,2)] = poly2cw(bp(:,1),bp(:,2));
% % % % % % %     bp = flipud(bp);
% % % % % % %     bp(end+1,:) = bp(1,:);
% % % % % % %     assert(ispolycw(bp(:,1),bp(:,2)));
% % % % % %     
% % % % % %     % cuboid bottom
% % % % % %     for j = 1:num_cuboid
% % % % % %         cp = cat(1, cuboidhyp(j).junc3([4 3 7 8]).pt);
% % % % % % %         assert(ispolycw(cp(:,1),cp(:,2)));
% % % % % %         [cp(:,1) cp(:,2)] = poly2cw(cp(:,1),cp(:,2));
% % % % % % %         cp = flipud(cp);
% % % % % % %         cp(end+1,:) = cp(1,:);
% % % % % % %         assert(ispolycw(cp(:,1),cp(:,2)));
% % % % % % 
% % % % % %         [x y] = polybool('subtraction', cp(:,1),cp(:,2), bp(:,1),bp(:,2));
% % % % % %         if ~isempty(x)
% % % % % %         [x y] = polybool('intersection', x, y, imgpoly(:,1), imgpoly(:,2));
% % % % % %         end
% % % % % %         
% % % % % %         if isempty(x)
% % % % % % %         if length(x)/size(cp,1) > 0.05
% % % % % %             compat(i,j) = true;
% % % % % %         else
% % % % % %             compat(i,j) = false;
% % % % % %         end
% % % % % %         
% % % % % % %         cuboid_bottom = poly2mask(pt(:,1), pt(:,2), imgsize(1), imgsize(2));
% % % % % % 
% % % % % % %         if any( cuboid_bottom(:)==1 & building_bottom(:)==0 )
% % % % % % %             compat(i,j) = false;
% % % % % % %         else
% % % % % % %             compat(i,j) = true;
% % % % % % %         end
% % % % % %     end
% % % % % % end

%%
function p = building_floor_poly(buildinghyp, imgsize)
% for i = 1:length(buildinghyp.corner)
%     assert(strcmp(buildinghyp.corner(i).type, 'minusorplus'));
% end

num_walls = length(buildinghyp.box);

num_vertex = num_walls+1;

if buildinghyp.box(1).p2(1) == 1 % if 1st wall ends at image left edge
    num_vertex = num_vertex+1;
end
if buildinghyp.box(end).p4(1) == imgsize(2) % if last wall ends at image right edge
    num_vertex = num_vertex+1;
end

p = zeros(num_vertex, 2);

for ivert = 1:num_walls
    p(ivert,:) = buildinghyp.box(ivert).p2;
end
ivert = ivert+1;
p(ivert,:) = buildinghyp.box(end).p4;

if buildinghyp.box(end).p4(1) == imgsize(2) % if last wall ends at image right edge
    ivert = ivert+1;
    p(ivert,:) = [imgsize(2) imgsize(1)];
end
if buildinghyp.box(1).p2(1) == 1 % if 1st wall ends at image left edge
    ivert = ivert+1;
    p(ivert,:) = [1 imgsize(1)];
end

