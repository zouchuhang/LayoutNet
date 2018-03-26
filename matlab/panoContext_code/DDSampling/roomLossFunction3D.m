function [ cost,gnd_objects,hyp_objects ] = roomLossFunction3D( gnd_objects, hyp_objects, H2G_R )
%ROOMLOSSFUNCTION Summary of this function goes here
%   Detailed explanation goes here

%% define cost
gnd_loss_cost = 2;
hyp_false_cost = 1;

%% convert hypothesis point
for i = 1:length(hyp_objects)
    hyp_objects(i).align = rotatePoint(hyp_objects(i).align, H2G_R);
end

%% greedy matching
test_type = [1:12 29];
test_type_num = length(test_type);
cost = 0;
hyp_objects(1).matchID = [];
hyp_objects(1).cost = [];

gnd_objects(1).matchID = [];
gnd_objects(1).cost = [];

for i = 1:test_type_num
    tid = test_type(i);
    
    hyp_rep_valid = [hyp_objects.objtype]==tid;
    hyp_rep_num = sum(hyp_rep_valid);
    hyp_rep = hyp_objects(hyp_rep_valid);
    hyp_rep_id = reshape(find(hyp_rep_valid),1,[]);
    
    gnd_rep_valid = [gnd_objects.objtype]==tid;
    gnd_rep_num = sum(gnd_rep_valid);
    gnd_rep = gnd_objects(gnd_rep_valid);
    gnd_rep_id = reshape(find(gnd_rep_valid),1,[]);
    
    if hyp_rep_num==0
        cost = cost + gnd_rep_num * gnd_loss_cost;
        for k = gnd_rep_id
            gnd_objects(k).cost = gnd_loss_cost;
            gnd_objects(k).matchID = 0;
        end
        
    elseif gnd_rep_num==0
        cost = cost + hyp_rep_num * hyp_false_cost;
        for k = hyp_rep_id
            hyp_objects(k).cost = hyp_false_cost;
            hyp_objects(k).matchID = 0;
        end
        
    else
        hyp_gnd_dist = zeros(hyp_rep_num, gnd_rep_num);
        hyp_gnd_maxs = zeros(hyp_rep_num, gnd_rep_num);
        for m = 1:hyp_rep_num
            for n = 1:gnd_rep_num
                
                hyp_gnd_dist(m,n) = verticesDistance(hyp_rep(m).align, gnd_rep(n).align);
                hyp_gnd_maxs(m,n) = 1;
%                 d = volumeIntersection(hyp_rep(m).point3D, gnd_rep(n).point3D);
%                 if d>=0
%                     hyp_gnd_dist(m,n) = 2*d;
%                     hyp_gnd_maxs(m,n) = 2;
%                 else
%                     hyp_gnd_dist(m,n) = verticesDistance(hyp_rep(m).point3D, gnd_rep(n).point3D);
%                     hyp_gnd_maxs(m,n) = 1;
%                 end
                
%                 [v1, v2, vi] = obbInteVolume(hyp_rep(m).point3D, gnd_rep(n).point3D);
%                 hyp_gnd_dist(m,n) = vi/(v1+v2-vi);
%                 p1 = hyp_rep(m).point3D;
%                 p2 = gnd_rep(n).point3D;
%                 meansize1 = (norm(p1(1,:)-p1(2,:)) + norm(p1(1,:)-p1(3,:)) + norm(p1(1,:)-p1(5,:)))/3;
%                 meansize2 = (norm(p2(1,:)-p2(2,:)) + norm(p2(1,:)-p2(3,:)) + norm(p2(1,:)-p2(5,:)))/3;
%                 meandist = mean(sqrt(sum((hyp_rep(m).point3D-gnd_rep(n).point3D).^2,2))) ...
%                     /sqrt(meansize1*meansize2);                
%                 hyp_gnd_dist(m,n) = sigmoidFunc(meandist, 0.45, 4); %exp(-2*meandist);%
                
%                 if hyp_gnd_dist(m,n)>0.5
%                     fprintf('type: %d, hypid: %d, gndid: %d, score: %f\n', tid, hyp_rep_id(m), gnd_rep_id(n), hyp_gnd_dist(m,n));
%                 end
                
                
            end
        end
        
        sel_hyp = false(hyp_rep_num,1);
        sel_gnd = false(gnd_rep_num,1);
        for k = 1:min( hyp_rep_num, gnd_rep_num)
            [r,c] = find(hyp_gnd_dist==max(hyp_gnd_dist(:)));
            r = r(1); c = c(1);
            cost = cost + (hyp_gnd_maxs(r,c)-hyp_gnd_dist(r,c));
            hyp_objects(hyp_rep_id(r)).cost = (hyp_gnd_maxs(r,c)-hyp_gnd_dist(r,c))/2;
            hyp_objects(hyp_rep_id(r)).matchID = gnd_rep_id(c);
            gnd_objects(gnd_rep_id(c)).cost = (hyp_gnd_maxs(r,c)-hyp_gnd_dist(r,c))/2;
            gnd_objects(gnd_rep_id(c)).matchID = hyp_rep_id(r);
            
            sel_hyp(r) = true;
            sel_gnd(c) = true;
            hyp_gnd_dist(:,c) = -1;
            hyp_gnd_dist(r,:) = -1;
        end
        
        if hyp_rep_num>gnd_rep_num
            cost = cost + (hyp_rep_num-gnd_rep_num)*hyp_false_cost;
        else
            cost = cost + (gnd_rep_num-hyp_rep_num)*gnd_loss_cost;
        end
        rr = find(~sel_hyp);
        rc = find(~sel_gnd);
        for k = 1:length(rr)
            hyp_objects(hyp_rep_id(rr(k))).cost = hyp_false_cost;
            hyp_objects(hyp_rep_id(rr(k))).matchID = 0;
        end
        for k = 1:length(rc)
            gnd_objects(gnd_rep_id(rc(k))).cost = gnd_loss_cost;
            gnd_objects(gnd_rep_id(rc(k))).matchID = 0;
        end
    end
    
end

cost = cost/length(gnd_objects);
end

function d = verticesDistance(p1, p2)
% p1 = hyp_rep(m).point3D;
% p2 = gnd_rep(n).point3D;
meansize1 = (norm(p1(1,:)-p1(2,:)) + norm(p1(1,:)-p1(4,:)) + norm(p1(1,:)-p1(5,:)))/3;
meansize2 = (norm(p2(1,:)-p2(2,:)) + norm(p2(1,:)-p2(4,:)) + norm(p2(1,:)-p2(5,:)))/3;
meandist = mean(sqrt(sum((p1-p2).^2,2))) ...
    /sqrt(meansize1*meansize2);                
d = sigmoidFunc(meandist, 0.45, 4); %exp(-2*meandist);%
end

% function d = volumeIntersection(p1, p2)
% [v1, v2, vi] = obbInteVolume(p1, p2);
% if v1==-1
%     d = -1;
% else
%     d = min(v1,v2)/vi;
% end
% end

function p = sigmoidFunc(v, m, a)
p = 1./(1+exp(a/m*(v-m)));
end

% gnd_objects_id = zeros(length(gnd_objects),1);
% gnd_valid = true(length(gnd_objects),1);
% gnd_xyz = zeros(length(gnd_objects),6);
% gnd_point = cell(length(gnd_objects),1);
% for i = 1:length(gnd_objects)
%     gnd_objects_id(i) = get_object_type({gnd_objects(i).name});
%     x = gnd_objects(i).x_w;
%     ad = abs(round(x(7)/(pi/2))*(pi/2) - x(7));
%     if ad>0.1
%         gnd_valid(i) = false;
%     end
%     out_points_w = gnd_objects(i).out_points_w;
%     gnd_xyz(i,:) = [min(out_points_w, [], 1) max(out_points_w, [], 1)];
%     gnd_point{i} = anno2point(out_points_w, annorule(rule_type==gnd_objects(i).type));
% end
% 
% 
% hyp_objects_id = zeros(length(hyp_objects),1);
% hyp_xyz = zeros(length(hyp_objects),6);
% hyp_point = zeros(8,3,length(hyp_objects));
% for i = 1:length(hyp_objects)
%     hyp_objects_id(i) = get_object_type({hyp_objects(i).name});
%     out_points_w = hyp_objects(i).out_points_w;
%     hyp_xyz(i,:) = [min(out_points_w, [], 1) max(out_points_w, [], 1)];
% end
% hyp_valid = true(length(hyp_objects),1);
% 
% subset_category = [1:12 29];
% XYZ2POINTLIST = [1 2 3; 4 2 3; 1 5 3; 4 5 3; ...
%                  1 2 6; 4 2 6; 1 5 6; 4 5 6];
% 
%              
% %% check precision
% categorycost = inf*ones(length(subset_category),1);
% categorynumb = zeros(length(subset_category),2);
% hyp_hit_id = zeros(length(hyp_objects),2);
% for k = 1:length(subset_category)
%     oid = subset_category(k);
%     gnd_ids = find(gnd_valid & gnd_objects_id==oid);   
%     hyp_ids = find(hyp_objects_id==oid);
%     num_gnd = length(gnd_ids);
%     num_hyp = length(hyp_ids);
%     categorynumb(k,:) = [num_hyp num_gnd];
%     if num_hyp == 0 && num_gnd~=0;
%         categorycost(k) = 860530;
%         continue;
%     elseif num_gnd == 0 && num_hyp~=0
%         categorycost(k) = 860530;
%         continue;
%     elseif num_gnd == 0 && num_hyp == 0
%         categorycost(k) = 0;
%         continue;
%     end
%     
%     distance = zeros(length(hyp_ids), length(gnd_ids));
%     gnd_pl = zeros(length(gnd_ids)*8,3);
%     for i = 1:length(gnd_ids)
%         xyz = gnd_xyz(gnd_ids(i),:);
%         gnd_pl(i*8-7:i*8,:) = xyz(XYZ2POINTLIST);
%     end
%     for i = 1:length(hyp_ids)
%         xyz = hyp_xyz(hyp_ids(i),:);
%         hyp_pl = xyz(XYZ2POINTLIST);
%         dist = sqrt(sum((repmat(hyp_pl, length(gnd_ids), 1) - gnd_pl).^2,2));
%         dist = sum(reshape(dist, 8, []),1);
%         distance(i,:) = dist;
%     end
% 
%     assigned_id = zeros(num_hyp, 1);
%     assigned_gnd = false(num_gnd,1);
%     assigned = false(num_hyp, 1);
%     assigned_dist = inf*ones(num_hyp, 1);
%     while ~all(assigned) && ~all(assigned_gnd)
%         [B,I] = min(distance, [], 2);
%         [C,J] = min(B);
%         assigned(J) = true;
%         assigned_gnd(I(J)) = true;
%         assigned_id(J) = I(J);
%         assigned_dist(J) = C;
%         hyp_hit_id(hyp_ids(J),:) = [gnd_ids(I(J)) C];
%         
%         distance(J,:) = inf;
%         distance(:,I(J)) = inf;
%     end
%     
%     max_penalty = 2*max(assigned_dist(~isinf(assigned_dist)));
%     assigned_dist(isinf(assigned_dist)) = max_penalty;
%     miss_gnd_num = max(0, num_gnd-num_hyp);
% %     try
%     categorycost(k) = mean([assigned_dist; max_penalty*ones(miss_gnd_num, 1)]);
% %     catch
% %         fprintf('');
% %     end
% end
% 
% cost = mean(categorycost(~isinf(categorycost)));
% 
% end
% 
% function points = anno2point(out_points, annorule)
% inv_vertex = annorule.inv_vertex;
% for i = 1:size(inv_vertex,1)
%     out_points(end+1,:) = out_points(inv_vertex(i,1),:) ...
%                         + out_points(inv_vertex(i,2),:) ...
%                         - out_points(inv_vertex(i,3),:);
% end   
% points = out_points;
% end

