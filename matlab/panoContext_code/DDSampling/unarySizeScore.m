function [ unaryscore, unaryanglesid ] = unarySizeScore( room_xyz, obj_xyz, typenum, unary_size)
%UNARYSIZESCORE Summary of this function goes here
%   Detailed explanation goes here
if isempty(obj_xyz)
    % a special test 
    if (room_xyz(5)-room_xyz(2))>(room_xyz(4)-room_xyz(1))
        data.testangle = 1;
    else
        data.testangle = 2;
    end
    room_unary_size = unary_size{29}; % make this as the room id
    swap = room_unary_size(:,1)>room_unary_size(:,2);
    t = room_unary_size(swap,2);
    room_unary_size(swap,2) = room_unary_size(swap,1);
    room_unary_size(swap,1) = t;

    a = mean(room_unary_size,1);
    p = sqrt(sum(a.^2));
    unaryscore = unaryHyperTest( room_xyz, room_unary_size, 1, data )/p;
    
else
    averagesize = zeros(typenum,3);
    for i = 1:typenum
        averagesize(i,:) = mean(unary_size{i},1);
    end
    trait = sqrt(sum(averagesize.^2,2));
    % unary_para  = trait(:,1)*15/trait(1,1);

    num_object = size(obj_xyz,1);
    unaryscore = zeros(num_object, typenum);
    unaryanglesid = zeros(num_object, 1);
    for oid = 1:num_object
        xyz = obj_xyz(oid,:);
        dist_diff = xyz - repmat(room_xyz, size(xyz,1), 1);
        [~, id] = min(abs(dist_diff(:,[1 2 4 5])),[],2);   
        unaryanglesid(oid) = id;
        data.roomrect = room_xyz;
        data.testangle = id;
        for tid = 1:typenum
            u = unaryHyperTest( xyz, unary_size{tid}, 1, data );
            unaryscore(oid,tid) = u/trait(tid); %sigmoidFunc(u, unary_para(tid), 5);
        end
    end
    
end

end

% function p = sigmoidFunc(v, m, a)
% p = 1./(1+exp(a/m*(v-m)));
% end
% 
% function p = sigmoidFunc_flex(v, a, b)
% p = 1./(1+exp(a*v+b));
% end