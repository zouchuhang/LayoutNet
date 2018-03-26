% folderName = 'annotation/dataset/bedroom/'; % buffer dataset
% load([folderName 'ANNO_ALL.mat']);
% load([folderName 'IMGLIST.mat']);
% data_valid = true(length(ANNO_ALL), 1);
% for i = 1:length(ANNO_ALL)
% %     fprintf('%d\n',i);
%     if isempty(ANNO_ALL(i).ANNO3D) || ~ANNO_ALL(i).ANNO3D.b_singleroom
%         data_valid(i) = false;
%     else
%         data_valid(i) = true;
%     end
% end
% VALID_SAMPLE_ID = find(data_valid);
% valid_sample_number = length(VALID_SAMPLE_ID);
% 
% %%
% total_fea_dim = comp_object_3d_feature();
% pairwise = cell(29,29);
% pairwise_norm = cell(29,29);
% pairwise_wall = cell(29,29);
% pairwise_dataid = cell(29,29);
% unary_size = cell(29,1);
% unary_size_norm = cell(29,1);
% % unary_posi = cell(29,1);
% % unary__norm = cell(29,1);
% unary_dataid = cell(29,1);
% 
% 
% for j = 1:valid_sample_number
%     fprintf('%d/%d\n', j, valid_sample_number);
%     did = VALID_SAMPLE_ID(j);
%     if isempty(ANNO_ALL(did).ANNO3D) || ~ANNO_ALL(did).ANNO3D.b_singleroom
%         continue;
%     end
%     objects3D = ANNO_ALL(did).ANNO3D.objects3D;
%     contain3Dcoords = ANNO_ALL(did).ANNO3D.contain3Dcoords;
%     room_id = contain3Dcoords(1);
%     room_obj = objects3D(room_id);
%     
%     
% %     unary_size{29,1} = [unary_size{29,1}; room_size];
% %     unary_size_norm{29,1} = [unary_size_norm{29,1}; [1 1 1]];
% %     unary_dataid{29,1} = [unary_dataid{29,1}; [did room_id]];
%     
%     for i = 1:length(contain3Dcoords)
%         aid = contain3Dcoords(i);
%         [objects angle] = alignObjectViaWall(objects3D, room_obj, aid, contain3Dcoords);
%         
%         room_xyz = [min(objects(room_id).out_points_w,[],1) max(objects(room_id).out_points_w,[],1)];
%         room_size = room_xyz(4:6)-room_xyz(1:3);
%         
%         lhs_points = objects(aid).out_points_w;
%         lhs_xyz = [min(lhs_points,[],1) max(lhs_points,[],1)];
%         lhs_size = (lhs_xyz(4:6)-lhs_xyz(1:3))+0.000001;
%         lhs_center = (lhs_xyz(4:6) + lhs_xyz(1:3))/2;
%         lhs_wall = [lhs_xyz(1) (lhs_xyz([2 3])+lhs_xyz([5 6]))/2];
%         lhsID = get_object_type({objects(aid).name});
%         
%         unary_size{lhsID,1} = [unary_size{lhsID};lhs_size];
%         unary_size_norm{lhsID,1} = [unary_size_norm{lhsID,1}; lhs_size./room_size];
%         unary_dataid{lhsID,1} = [unary_dataid{lhsID,1}; [did aid]];
%         
%         for k = setdiff(contain3Dcoords, aid)
%             rhsID = get_object_type({objects(k).name});
%             rhs_points = objects(k).out_points_w;
%             rhs_xyz = [min(rhs_points,[],1) max(rhs_points,[],1)];
%             rhs_center = (rhs_xyz(4:6) + rhs_xyz(1:3))/2;
%             rhs_wall = [rhs_xyz(1) (rhs_xyz([2 3])+rhs_xyz([5 6]))/2];
%             
%             pairwise{lhsID, rhsID} = [pairwise{lhsID,rhsID}; (rhs_center-lhs_center)];
%             pairwise_norm{lhsID,rhsID} = ...
%                 [pairwise_norm{lhsID,rhsID}; (rhs_center-lhs_center)./lhs_size];
%             pairwise_wall{lhsID, rhsID} = [pairwise_wall{lhsID, rhsID}; rhs_wall-lhs_wall];
%             pairwise_dataid{lhsID, rhsID} = [pairwise_dataid{lhsID, rhsID}; [did aid k]];
%         end
%         
%     end
% end
% 
% save('object_classifier/objectpairwise.mat','pairwise','pairwise_norm','pairwise_wall','pairwise_dataid');
% save('object_classifier/objectunary.mat','unary_size','unary_size_norm','unary_dataid');
%%
load('object_classifier/objectpairwise.mat');
[m,n] = size(pairwise);

for i = 1:m
    for j = 1:n
%         figure(1);  clf;
%         pointnum = size(pairwise{i,j},1);
%         points = pairwise{i,j};
%         subplot(1,2,1); 
%         for k = 1:pointnum
%             plot3(points(k,1), points(k,2), points(k,3), ...
%                 '--rs',...
%                 'MarkerFaceColor','g',...
%                 'MarkerSize',2); hold on;
%         end
%         view(2);
        if isempty(pairwise{i,j})
            continue;
        end

        pts = pairwise{i,j};
        pts_add = pts;
        pts_add(:,2) = pts_add(:,2)*-1;
        pairwise{i,j} = [pts;pts_add];
        
        
%         pointnum = size(pairwise{i,j},1);
%         points = pairwise{i,j};
%         subplot(1,2,2); 
%         for k = 1:pointnum
%             plot3(points(k,1), points(k,2), points(k,3), ...
%                 '--rs',...
%                 'MarkerFaceColor','g',...
%                 'MarkerSize',2);  hold on;
%         end
%         view(2);
        
        pts = pairwise_norm{i,j};
        pts_add = pts;
        pts_add(:,2) = pts_add(:,2)*-1;
        pairwise_norm{i,j} = [pts;pts_add];
        
        pts = pairwise_wall{i,j};
        pts_add = pts;
        pts_add(:,2) = pts_add(:,2)*-1;
        pairwise_wall{i,j} = [pts;pts_add];
        
        pts = pairwise_dataid{i,j};
        pts = [pts zeros(size(pts,1),1)];
        pts_add = pts;
        pts_add(:,4) = 1;
        pairwise_dataid{i,j} = [pts;pts_add];
    end
end
save('object_classifier/objectpairwise_flip.mat','pairwise','pairwise_norm','pairwise_wall','pairwise_dataid');

%%

mkdir('pairwise_wall');
for lhs = [1:12 29]
    for rhs = [1:12 29]
        figure(1);  clf;
        % lhs = 9;
        lhs_name = get_object_type(lhs);
        % rhs = 2;
        rhs_name = get_object_type(rhs);
        points = pairwise{lhs,rhs};
        points_norm = pairwise_wall{lhs,rhs};
        pointnum = size(points,1);
        subplot(1,2,1); 
        for i = 1:pointnum
            plot3(points(i,1), points(i,2), points(i,3), ...
                '--rs',...
                'MarkerFaceColor','g',...
                'MarkerSize',2);  hold on;
        end
        title(sprintf('%s->%s', lhs_name{1}, rhs_name{1}));
        axis equal; axis([-200 800 -500 500 -500 500]); view(2);

        subplot(1,2,2); 
        for i = 1:pointnum
            plot3(points_norm(i,1), points_norm(i,2), points_norm(i,3), ...
                '--rs',...
                'MarkerFaceColor','g',...
                'MarkerSize',2);  hold on;
        end
        title(sprintf('%s->%s', lhs_name{1}, rhs_name{1}));
        axis equal; axis([-200 800 -500 500 -500 500]); view(2);
%         axis equal; axis([-2 10 -6 6 -6 6]); view(2); 

        print('-djpeg',sprintf('./pairwise_wall/%s_%s.jpg', lhs_name{1}, rhs_name{1}));
    end
end

%%
for lhs = 9
    for rhs = 1
        figure(2);  clf;
        % lhs = 9;
        lhs_name = get_object_type(lhs);
        % rhs = 2;
        rhs_name = get_object_type(rhs);
%         points = pairwise;
        points = rotate_pairwise{lhs,rhs,3};
        points_norm = rotate_pairwise{lhs,rhs,4};
        pointnum = size(points,1);
        subplot(1,2,1); 
        for i = 1:pointnum
            plot3(points(i,1), points(i,2), points(i,3), ...
                '--rs',...
                'MarkerFaceColor','g',...
                'MarkerSize',2);  hold on;
        end
        title(sprintf('%s->%s', lhs_name{1}, rhs_name{1}));
        axis equal;  view(2);

        subplot(1,2,2); 
        for i = 1:pointnum
            plot3(points_norm(i,1), points_norm(i,2), points_norm(i,3), ...
                '--rs',...
                'MarkerFaceColor','g',...
                'MarkerSize',2);  hold on;
        end
        title(sprintf('%s->%s', lhs_name{1}, rhs_name{1}));
        axis equal;  view(2); 

%         print('-djpeg',sprintf('./pairwise/%s_%s.jpg', lhs_name{1}, rhs_name{1}));
    end
end

%%
points = rotate_pairwise_wall{i,tid,type_hyps(i).anglesid(j)};
figure(3);  clf;
pointnum = size(points,1);

for iii = 1:pointnum
    plot3(points(iii,1), points(iii,2), points(iii,3), ...
        '--rs',...
        'MarkerFaceColor','g',...
        'MarkerSize',2);  hold on;
end
axis equal;  view(2);


