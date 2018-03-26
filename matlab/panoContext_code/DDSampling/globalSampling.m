function [ ALL_HYPS ] = globalSampling( hyps_obj, number, typenum )
%GLOBALSAMPLING Summary of this function goes here
%   Detailed explanation goes here
% room_xyz = [min(room_obj.out_points_w,[],1) max(room_obj.out_points_w,[],1)];
% typenum = 12;
global config

%% rotate pairwise constraints
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

dpmvalid = true(length(rectscore),typenum);
% dpmscore = vertcat(OBJFEA.dpmscore);
% dpmvalid(:,10) = dpmscore(:,1)>-0.9;
% dpmvalid(:,12) = dpmscore(:,2)>-0.95;

imgvalid = rectscore.*segmscore>0.1;
semvalid = unaryscore(:,1:typenum).*randomforest(:,1:typenum)>0.3;

totalvalid = dpmvalid & semvalid & repmat(imgvalid, 1, typenum);
totalid = cell(typenum,1);
for i = 1:typenum
    totalid{i} = find(totalvalid(:,i));
end

% rectscore = sqrt(rectscore);
% scr10 =  sigmoidFunc_flex(dpmscore(:,1), -50, -0.9);
% scr12 =  sigmoidFunc_flex(dpmscore(:,2), -15, -0.9);
% unaryscore(:,10) = unaryscore(:,10).*scr10;
% unaryscore(:,12) = unaryscore(:,12).*scr12;

%% decide object sequence and number
numberbias = exp(0.02*(-1:-1:-10));
load(config.numberName);
ALL_TYPERANK = cell(number,1);
for setid = 1:number
    sample_score = zeros(1,0);
    sample_ids = zeros(1,0);
    num_of_type = zeros(1,typenum);
    for oid = 1:typenum
        num_of_type(oid) = randp( prob(oid,1:10).*numberbias, 1) - 1;
        likelihood = rectscore.*segmscore.*unaryscore(:,oid).*randomforest(:,oid);
        sample_score = [sample_score likelihood(randp(likelihood, 1, num_of_type(oid)))'];
        sample_ids = [sample_ids oid*ones(1, num_of_type(oid))];
    end
    [~, I] = sort(sample_score, 'descend');
    ALL_TYPERANK{setid} = sample_ids(I);
end

%% start sampling
occ_thres = 0.05;
par_thres = [50 0.3 50];
ori_sel_hyps = repmat(struct('fixed',false,'selID',[]), typenum, 1);

ALL_HYPS = repmat(struct('sel_hyps',[]), number, 1);
for count = 1:number % parfor
    fprintf('Whole Room Hypothesis: %d/%d\n', count, number);
    typerank = ALL_TYPERANK{count};   
    sel_hyps = ori_sel_hyps;
   
    degree = 2*ones(1,length(typerank));
    
    for samid = 1:length(typerank)
        tid = typerank(samid);
        unary_semantic = unaryscore(totalvalid(:,tid),tid).*randomforest(totalvalid(:,tid),tid); % type_hyps(tid).unary_score; %type_hyps(tid).h_scores(:, tid); %
        unary_imageevi = rectscore(totalvalid(:,tid)).*segmscore(totalvalid(:,tid)); % type_hyps(tid).image_score; %ones(length(unary_semantic),1); %
        pairwise_abso = ones(length(unary_semantic),1);
        pairwise_norm = ones(length(unary_semantic),1);
        pairwise_wall = ones(length(unary_semantic),1);  
        pairwise_occl = ones(length(unary_semantic),1);
        
        % check fixed objects and remove invalid hypothesis
        xyz = hyps_obj.obj_xyz(totalvalid(:,tid),:);
        for i = 1:typenum            
            if sel_hyps(i).fixed
                for j = sel_hyps(i).selID
                    if ~((tid==7 && i==10) || (i==7 && tid==10))
                        o = occlusionTest(hyps_obj.obj_xyz(j,:), xyz);
                    else
                        o = zeros(size(xyz,1),1);
                    end
%                     locdata = data;
                    locdata.seedangle = hyps_obj.anglesid(j);
                    locdata.testangle = hyps_obj.anglesid(totalvalid(:,tid));
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
        sub_prob = unary_semantic.*unary_imageevi.*pairwise_abso.*pairwise_wall.*pairwise_occl;
        if sum(sub_prob.^degree(samid))==0
            continue;
        end
        
        try
            sub_sel = totalid{tid}(randp(sub_prob.^degree(samid), 1));
        catch
            continue;
        end
        
        sel_hyps(tid).fixed = true;
        sel_hyps(tid).selID = [sel_hyps(tid).selID sub_sel];
        
    end
    
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
