function [ ALL_HYPS ] = addObject( hyps_obj, INP_HYPS, TPS, whitelist, scenescore )
%GLOBALSAMPLING Summary of this function goes here
%   INP_HYPS: list of hypothesis waiting for adding
%   TPS: type of object added for each hypothesis
%   whitelist: each column assign hard whitelist for each hypothesis
% room_xyz = [min(room_obj.out_points_w,[],1) max(room_obj.out_points_w,[],1)];
typenum = size(scenescore,2);

%% rotate pairwise constraints
global config
% load('./object_classifier/objectpairwise.mat');
load(config.pairwiseName);

angleset = [0 -pi/2 pi +pi/2];
rotate_pairwise = cell(size(pairwise,1), size(pairwise,2), 4);
rotate_pairwise_norm = cell(size(pairwise_norm,1), size(pairwise_norm,2), 4);
rotate_pairwise_wall = cell(size(pairwise_wall,1), size(pairwise_wall,2), 4);

for a = 1:4
    for i = 1:size(pairwise,1)
        for j = 1:size(pairwise,2)
            if ~isempty(pairwise{i,j})
                p = pairwise{i,j};
                r = p;
                r(:,1) = p(:,1).*cos(angleset(a)) + p(:,2).*sin(angleset(a));
                r(:,2) = -p(:,1).*sin(angleset(a)) + p(:,2).*cos(angleset(a));
                rotate_pairwise{i,j,a} = r;
                
                p = pairwise_norm{i,j};
                r = p;
                r(:,1) = p(:,1).*cos(angleset(a)) + p(:,2).*sin(angleset(a));
                r(:,2) = -p(:,1).*sin(angleset(a)) + p(:,2).*cos(angleset(a));
                rotate_pairwise_norm{i,j,a} = r;
                
                p = pairwise_wall{i,j};
                r = p;
                r(:,1) = p(:,1).*cos(angleset(a)) + p(:,2).*sin(angleset(a));
                r(:,2) = -p(:,1).*sin(angleset(a)) + p(:,2).*cos(angleset(a));
                rotate_pairwise_wall{i,j,a} = r;
            end
        end
    end
end

%% convert every score to 0~1 by sigmoid function
OBJFEA = hyps_obj.objfea;
unaryscore = vertcat(OBJFEA.unaryscore);
rectscore = vertcat(OBJFEA.rectscore);
segmscore = vertcat(OBJFEA.segmscore);
randomforest = vertcat(OBJFEA.randomforest);

unaryscore = sigmoidFunc_flex(unaryscore, 33.67, 0.15);
rectscore = sigmoidFunc_flex(rectscore, -4, 0);
segmscore = ones(length(segmscore),1);
% segmscore = sigmoidFunc_flex(segmscore, -4, 0.5);
randomforest = sigmoidFunc_flex(randomforest, -4, 0.5);

valid1 = rectscore.*segmscore>0.1;
valid2 = unaryscore(:,1:typenum).*randomforest(:,1:typenum)>0.3;

% valid1 = rectscore>0.1;
% valid2 = randomforest(:,1:typenum)>0.1;

priorvalid = repmat(valid1,1,typenum) & valid2;
% priorid = cell(typenum,1);
% for i = 1:12
%     priorid{i} = find(priorvalid(:,i));
% end

%% start sampling
occ_thres = 0.05;
par_thres = [50 0.3 50];

ALL_HYPS = INP_HYPS;
number = length(INP_HYPS);
for count = 1:number % parfor
%     fprintf('Whole Room Hypothesis: %d/%d\n', count, number);
    
    sel_hyps = ALL_HYPS(count).sel_hyps;
    
    tid = TPS(count);
    totalvalid = priorvalid(:,tid) & whitelist(:,count);
    totalid = find(totalvalid);
        
    unary_semantic = unaryscore(totalvalid,tid).*randomforest(totalvalid,tid); 
    unary_imageevi = rectscore(totalvalid).*segmscore(totalvalid); 
    unary_sceneevi = scenescore(totalvalid,tid);
    pairwise_abso = ones(length(unary_semantic),1);
    pairwise_norm = ones(length(unary_semantic),1);
    pairwise_wall = ones(length(unary_semantic),1);  
    pairwise_occl = ones(length(unary_semantic),1);

    % check fixed objects and remove invalid hypothesis
    xyz = hyps_obj.obj_xyz(totalvalid,:);
    for i = 1:typenum            
        if sel_hyps(i).fixed
            for j = sel_hyps(i).selID
                if ~((tid==7 && i==10) || (i==7 && tid==10))
                    o = occlusionTest(hyps_obj.obj_xyz(j,:), xyz);
                else
                    o = zeros(size(xyz,1),1);
                end

                locdata.seedangle = hyps_obj.anglesid(j);
                locdata.testangle = hyps_obj.anglesid(totalvalid);
                p1 = pairwiseHyperTest( hyps_obj.obj_xyz(j,:), xyz, ...
                    rotate_pairwise{i,tid,hyps_obj.anglesid(j)}, ...
                    1, i, tid, locdata);
%                     p2 = pairwiseHyperTest( type_hyps(i).h_obj_xyz(j,:), xyz, ...
%                         rotate_pairwise_norm{i,tid,type_hyps(i).anglesid(j)}, ...
%                         2, i, tid, data);
                p3 = pairwiseHyperTest( hyps_obj.obj_xyz(j,:), xyz, ...
                    rotate_pairwise_wall{i,tid,hyps_obj.anglesid(j)}, ...
                    3, i, tid, locdata);

                pairwise_occl = pairwise_occl .* (o<occ_thres);
                pairwise_abso = pairwise_abso .* sigmoidFunc(p1, par_thres(1), 5);
%                     pairwise_norm = pairwise_norm .* sigmoidFunc(p2, 0.3, 5);
                pairwise_wall = pairwise_wall .* sigmoidFunc(p3, par_thres(3), 5);
%                     validID = validID & o<occ_thres & p1<par_thres(1) & p2<par_thres(2) & p3<par_thres(3);
            end
        end
    end

    % randomly pick one              
    sub_prob = unary_sceneevi.*unary_semantic.*unary_imageevi.*pairwise_abso.*pairwise_wall.*pairwise_occl;
    if sum(sub_prob)==0
        continue;
    end
    sub_sel = totalid(randp(sub_prob, 1));

    sel_hyps(tid).fixed = true;
    sel_hyps(tid).selID = [sel_hyps(tid).selID sub_sel];
        
    ALL_HYPS(count).sel_hyps = sel_hyps;
end

end

%%
function p = sigmoidFunc(v, m, a)
p = 1./(1+exp(a/m*(v-m)));
end

function p = sigmoidFunc_flex(v, a, b)
% p = 1./(1+exp(a*v+b));
p = 1./(1+exp(a*(v-b)));
end

