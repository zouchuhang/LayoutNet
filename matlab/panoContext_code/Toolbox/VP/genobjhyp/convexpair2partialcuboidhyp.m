function partcubhyp = convexpair2partialcuboidhyp(convexpair, boundingquad)

[r c v] = find(convexpair);

partcubhyp = [];
for i = 1:length(r)
    bq1 = boundingquad(r(i));
    bq2 = boundingquad(c(i));
    conftype = v(i);
    partcubhyp = [partcubhyp ...
        convexpair2partialcuboidhyp_single(bq1, bq2, conftype)];
end

%%
function pch = convexpair2partialcuboidhyp_single(bq1, bq2, conftype)

junc3dummy = struct('type',[],'pt',[],'regionid',[]);
c = 0;

% 1 if top & front & region3or4
if     conftype == 1 
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(4).pt = bq2.pt(4,:);
    pch(c).junc3(5).pt = bq1.pt(1,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(3).pt = bq2.pt(3,:);
    pch(c).junc3(5).pt = bq1.pt(1,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);

    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(3).pt = bq2.pt(3,:);
    pch(c).junc3(4).pt = bq2.pt(4,:);
    pch(c).junc3(5).pt = bq1.pt(1,:);

    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(3).pt = bq2.pt(3,:);
    pch(c).junc3(4).pt = bq2.pt(4,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);

    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(1).pt = bq1.pt(3,:);
    pch(c).junc3(4).pt = bq2.pt(4,:);
    pch(c).junc3(5).pt = bq1.pt(1,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(2).pt = bq1.pt(4,:);
    pch(c).junc3(3).pt = bq2.pt(3,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(1).pt = bq2.pt(1,:);
    pch(c).junc3(3).pt = bq2.pt(3,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(2).pt = bq2.pt(2,:);
    pch(c).junc3(4).pt = bq2.pt(4,:);
    pch(c).junc3(5).pt = bq1.pt(1,:);
    
% 2 if top & left & region4
elseif conftype == 2
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(1).pt = bq1.pt(3,:);
    pch(c).junc3(2).pt = bq1.pt(4,:);
    pch(c).junc3(7).pt = bq2.pt(3,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(5).pt = bq1.pt(1,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);
    pch(c).junc3(3).pt = bq2.pt(4,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(2).pt = bq1.pt(4,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);
    pch(c).junc3(3).pt = bq2.pt(4,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(2).pt = bq1.pt(4,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);
    pch(c).junc3(7).pt = bq2.pt(3,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(1).pt = bq2.pt(2,:);
    pch(c).junc3(3).pt = bq2.pt(4,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(2).pt = bq1.pt(4,:);
    pch(c).junc3(5).pt = bq2.pt(1,:);
    pch(c).junc3(7).pt = bq2.pt(3,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(3).pt = bq2.pt(4,:);
    pch(c).junc3(7).pt = bq2.pt(3,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(3).pt = bq2.pt(4,:);
    pch(c).junc3(7).pt = bq2.pt(3,:);
    pch(c).junc3(2).pt = bq1.pt(4,:);
    
% 3 if top & right & region3
elseif conftype == 3
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(1).pt = bq1.pt(3,:);
    pch(c).junc3(2).pt = bq1.pt(4,:);
    pch(c).junc3(8).pt = bq2.pt(4,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(5).pt = bq1.pt(1,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);
    pch(c).junc3(4).pt = bq2.pt(3,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(1).pt = bq1.pt(3,:);
    pch(c).junc3(5).pt = bq1.pt(1,:);
    pch(c).junc3(4).pt = bq2.pt(3,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(1).pt = bq1.pt(3,:);
    pch(c).junc3(5).pt = bq1.pt(1,:);
    pch(c).junc3(8).pt = bq2.pt(4,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(5).pt = bq1.pt(1,:);
    pch(c).junc3(2).pt = bq2.pt(1,:);
    pch(c).junc3(4).pt = bq2.pt(3,:);

    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(1).pt = bq1.pt(3,:);
    pch(c).junc3(6).pt = bq2.pt(2,:);
    pch(c).junc3(8).pt = bq2.pt(4,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(5).pt = bq1.pt(1,:);
    pch(c).junc3(4).pt = bq2.pt(3,:);
    pch(c).junc3(4).pt = bq2.pt(3,:);
    
% 4 if left & front & region4
elseif conftype == 4
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(1).pt = bq2.pt(1,:);
    pch(c).junc3(2).pt = bq2.pt(2,:);
    pch(c).junc3(7).pt = bq1.pt(3,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(3).pt = bq2.pt(3,:);
    pch(c).junc3(4).pt = bq2.pt(4,:);
    pch(c).junc3(5).pt = bq1.pt(1,:);

    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(2).pt = bq2.pt(2,:);
    pch(c).junc3(4).pt = bq2.pt(4,:);
    pch(c).junc3(5).pt = bq1.pt(1,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(2).pt = bq2.pt(2,:);
    pch(c).junc3(4).pt = bq2.pt(4,:);
    pch(c).junc3(7).pt = bq1.pt(3,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(1).pt = bq1.pt(2,:);
    pch(c).junc3(5).pt = bq1.pt(1,:);
    pch(c).junc3(4).pt = bq2.pt(4,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(3).pt = bq1.pt(4,:);
    pch(c).junc3(7).pt = bq1.pt(3,:);
    pch(c).junc3(2).pt = bq2.pt(2,:);

    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(2).pt = bq2.pt(2,:);
    pch(c).junc3(5).pt = bq1.pt(1,:);
    pch(c).junc3(7).pt = bq1.pt(3,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(4).pt = bq2.pt(4,:);
    pch(c).junc3(5).pt = bq1.pt(1,:);
    pch(c).junc3(7).pt = bq1.pt(3,:);
    
% 5 if right & front & region3
elseif conftype == 5
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(1).pt = bq2.pt(1,:);
    pch(c).junc3(2).pt = bq2.pt(2,:);
    pch(c).junc3(8).pt = bq1.pt(4,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(3).pt = bq2.pt(3,:);
    pch(c).junc3(4).pt = bq2.pt(4,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);

    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(1).pt = bq2.pt(1,:);
    pch(c).junc3(3).pt = bq2.pt(3,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(1).pt = bq2.pt(1,:);
    pch(c).junc3(3).pt = bq2.pt(3,:);
    pch(c).junc3(8).pt = bq1.pt(4,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(3).pt = bq2.pt(3,:);
    pch(c).junc3(2).pt = bq1.pt(1,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(1).pt = bq2.pt(1,:);
    pch(c).junc3(4).pt = bq1.pt(3,:);
    pch(c).junc3(8).pt = bq1.pt(4,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(1).pt = bq2.pt(1,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);
    pch(c).junc3(8).pt = bq1.pt(4,:);
    
    c = c+1;
    pch(c).junc3(8) = junc3dummy;
    pch(c).junc3(3).pt = bq2.pt(3,:);
    pch(c).junc3(6).pt = bq1.pt(2,:);
    pch(c).junc3(8).pt = bq1.pt(4,:);
end


    
%%
% function pch = convexpair2partialcuboidhyp_single_notminimumset(bq1, bq2, conftype)
% 
%     partcubhyp(i).junc3(8) = junc3dummy; % initialize
% 
%     % 1 if top & front & region3or4
%     if     conftype == 1 
%         partcubhyp(i).junc3(1).pt = mean([bq1.pt(3,:); ...
%                                           bq2.pt(1,:)],1);
%         partcubhyp(i).junc3(2).pt = mean([bq1.pt(4,:); ...
%                                           bq2.pt(2,:)],1);
%         partcubhyp(i).junc3(3).pt = bq2.pt(3,:);
%         partcubhyp(i).junc3(4).pt = bq2.pt(4,:);
%         partcubhyp(i).junc3(5).pt = bq1.pt(1,:);
%         partcubhyp(i).junc3(6).pt = bq1.pt(2,:);
%     % 2 if top & left & region4
%     elseif conftype == 2
%         partcubhyp(i).junc3(1).pt = mean([bq1.pt(3,:); ...
%                                          bq2.pt(2,:)],1);
%         partcubhyp(i).junc3(2).pt = bq1.pt(4,:);
%         partcubhyp(i).junc3(3).pt = bq2.pt(4,:);
%         partcubhyp(i).junc3(5).pt = mean([bq1.pt(1,:); ...
%                                          bq2.pt(1,:)], 1);
%         partcubhyp(i).junc3(6).pt = bq1.pt(2,:);
%         partcubhyp(i).junc3(7).pt = bq2.pt(3,:);
%     % 3 if top & right & region3
%     elseif conftype == 3
%         partcubhyp(i).junc3(1).pt = bq1.pt(3,:);
%         partcubhyp(i).junc3(2).pt = mean([bq1.pt(4,:); ...
%                                          bq2.pt(1,:)], 1);
%         partcubhyp(i).junc3(4).pt = bq2.pt(3,:);
%         partcubhyp(i).junc3(5).pt = bq1.pt(1,:);
%         partcubhyp(i).junc3(6).pt = mean([bq1.pt(2,:); ...
%                                          bq2.pt(2,:)], 1);
%         partcubhyp(i).junc3(8).pt = bq2.pt(4,:);
%     % 4 if left & front & region4
%     elseif conftype == 4
%         partcubhyp(i).junc3(1).pt = mean([bq1.pt(2,:); ...
%                                          bq2.pt(1,:)], 1);
%         partcubhyp(i).junc3(2).pt = bq2.pt(2,:);
%         partcubhyp(i).junc3(3).pt = mean([bq1.pt(4,:); ...
%                                          bq2.pt(3,:)], 1);
%         partcubhyp(i).junc3(4).pt = bq2.pt(4,:);
%         partcubhyp(i).junc3(5).pt = bq1.pt(1,:);
%         partcubhyp(i).junc3(7).pt = bq1.pt(3,:);
%     % 5 if right & front & region3
%     elseif conftype == 5
%         partcubhyp(i).junc3(1).pt = bq2.pt(1,:);
%         partcubhyp(i).junc3(2).pt = mean([bq1.pt(1,:); ...
%                                          bq2.pt(2,:)], 1);
%         partcubhyp(i).junc3(3).pt = bq2.pt(3,:);
%         partcubhyp(i).junc3(4).pt = mean([bq1.pt(3,:); ...
%                                          bq2.pt(4,:)], 1);
%         partcubhyp(i).junc3(6).pt = bq1.pt(2,:);
%         partcubhyp(i).junc3(8).pt = bq1.pt(4,:);
%     end
% 
