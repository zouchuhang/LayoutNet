function [ ALL_HYPS ] = globalSamplingWhitelist( hyps_obj, room_obj, number, typenum, whitelist, scenescore )
%GLOBALSAMPLING Summary of this function goes here
%   Detailed explanation goes here
% room_xyz = [min(room_obj.out_points_w,[],1) max(room_obj.out_points_w,[],1)];
% typenum = 12;

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

imgvalid = rectscore.*segmscore>0.1;
semvalid = unaryscore(:,1:typenum).*randomforest(:,1:typenum)>0.3;

totalvalid = semvalid & repmat(imgvalid, 1, typenum) & whitelist;
totalid = cell(typenum,1);
for i = 1:12
    totalid{i} = find(totalvalid(:,i));
end

%% decide object sequence and number
% load('./object_classifier/objectnumber.mat');
load(config.numberName);
ALL_TYPERANK = zeros(number, typenum);
ALL_TYPENUMB = zeros(number, typenum);
for setid = 1:number
    sample_score = zeros(1,typenum);
    num_of_type = zeros(1,typenum);
    for oid = 1:typenum
        if isempty(totalid{oid})
            sample_score(oid) = 0;
            num_of_type(oid) = 0;
        else
            likelihood = rectscore.*segmscore.*unaryscore(:,oid).*randomforest(:,oid);%.*scenescore(:,oid);
            s = sort(likelihood,'descend');
            s = s(1:min(length(s),10));
            numberprior = -1*ones(1,16);
            numberprior(1) = prob(oid,1)*(1-s(1));           
            for i = 1:length(s)  
                numberprior(i+1) = prob(oid, i+1)*mean(s(1:i));
            end
            
            numberprior(numberprior==-1) = [];
            num_of_type(oid) = randp(numberprior, 1) - 1;
            
%             s = type_hyps(oid).unary_score;%.*sigmoidFunc_flex(type_hyps(oid).h_bottomscore,-5,1.25);
            s = unaryscore(:,oid).*randomforest(:,oid);
            sample_score(oid) = s(randp(s,1));
        end
    end
    [~,typerank] = sort(sample_score,'descend');
    ALL_TYPERANK(setid,:) = typerank;
    ALL_TYPENUMB(setid,:) = num_of_type;
end

%% start sampling
occ_thres = 0.05;
par_thres = [50 0.3 50];
ori_sel_hyps = repmat(struct('fixed',false,'selID',[]), typenum, 1);

% ALL_HYPS = repmat(struct('objIDs',[],'objTPs',[]), number, 1); %cell(number,1);
ALL_HYPS = repmat(struct('sel_hyps',[]), number, 1);
for count = 1:number % parfor
%     count = count + 1;
%     fprintf('Whole Room Hypothesis: %d/%d\n', count, number);
    
    num_of_type = ALL_TYPENUMB(count,:);
    typerank = ALL_TYPERANK(count,:);
%     for tid = 1:typenum
%         sel_hyps(tid).fixed = false;
%         sel_hyps(tid).selID = [];
%     end
    sel_hyps = ori_sel_hyps;
    
    for tid = typerank
        if num_of_type(tid)==0
            continue;
        end
        
        unary_semantic = unaryscore(totalvalid(:,tid),tid).*randomforest(totalvalid(:,tid),tid); % type_hyps(tid).unary_score; %type_hyps(tid).h_scores(:, tid); %
        unary_imageevi = rectscore(totalvalid(:,tid)).*segmscore(totalvalid(:,tid)); % type_hyps(tid).image_score; %ones(length(unary_semantic),1); %
        unary_sceneevi = scenescore(totalvalid(:,tid),oid);
        pairwise_abso = ones(length(unary_semantic),1);
        pairwise_norm = ones(length(unary_semantic),1);
        pairwise_wall = ones(length(unary_semantic),1);  
        pairwise_occl = ones(length(unary_semantic),1);
        
        % check fixed objects and remove invalid hypothesis
%         xyz = type_hyps(tid).h_obj_xyz; % all possible IDs
        xyz = hyps_obj.obj_xyz(totalvalid(:,tid),:);
%         validID = true(size(xyz,1),1);
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
                    
                    if any(isnan(p1))
                        p1 = zeros(sum(totalvalid),1);
                        p3 = zeros(sum(totalvalid),1);
                    end
                    
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
        sub_sel = totalid{tid}(randp(sub_prob, 1));
        
        for i = 2:num_of_type(tid)          
            o = occlusionTest(hyps_obj.obj_xyz(sub_sel(i-1),:), xyz);
%             locdata = data;
            locdata.seedangle = hyps_obj.anglesid(sub_sel(i-1));
            locdata.testangle = hyps_obj.anglesid(totalvalid(:,tid));
            p1 = pairwiseHyperTest( hyps_obj.obj_xyz(sub_sel(i-1),:), xyz, ...
                rotate_pairwise{tid, tid, locdata.seedangle}, ...
                1, tid, tid, locdata);
%             p2 = pairwiseHyperTest( xyz(sub_sel(i-1),:), xyz, ...
%                 rotate_pairwise_norm{tid, tid, data.seedangle}, ...
%                 2, tid, tid, data);
            p3 = pairwiseHyperTest( hyps_obj.obj_xyz(sub_sel(i-1),:), xyz, ...
                rotate_pairwise_wall{tid, tid, locdata.seedangle}, ...
                3, tid, tid, locdata);

            pairwise_occl = pairwise_occl .* (o<occ_thres);
            pairwise_abso = pairwise_abso .* sigmoidFunc(p1, par_thres(1), 5);
%                     pairwise_norm = pairwise_norm .* sigmoidFunc(p2, 0.3, 5);
            pairwise_wall = pairwise_wall .* sigmoidFunc(p3, par_thres(3), 5);
            sub_prob = unary_semantic.*unary_imageevi.*pairwise_abso.*pairwise_wall.*pairwise_occl;
            if sum(sub_prob)==0
                break;
            else                
                sub_sel(1,end+1) = totalid{tid}(randp( sub_prob,1));
            end
        end
        
        sel_hyps(tid).fixed = true;
        sel_hyps(tid).selID = sub_sel;
        
    end
    
%     objIDs = zeros(0,1);
%     objTPs = zeros(0,1);
%     for i = 1:length(sel_hyps)
%         if sel_hyps(i).fixed
%             objIDs = [objIDs; sel_hyps(i).selID(:)];
%             objTPs = [objTPs; i*ones(length(sel_hyps(i).selID),1)];
%         end
%     end
%     ALL_HYPS(count).objIDs = objIDs;
%     ALL_HYPS(count).objTPs = objTPs;
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

