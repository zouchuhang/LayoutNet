function [ objects ] = rectangleCombo3Hypothesis( rects, views, room_points )
%RECTANGLECOMBO3HYPOTHESIS Summary of this function goes here
%   Detailed explanation goes here
normal = [-1 0 0; 1 0 0; 0 -1 0; 0 1 0; 0 0 -1; 0 0 1];
rects(views==0,:) = [];
rank_rects = zeros(size(rects));
views(views==0) = [];
for vid = 1:length(views)
    vp = normal(views(vid),:);
    xyz = reshape(rects(vid,:), 3, 4)';
    [out2D, valid, division] = projectPoint2SeparateView( xyz, vp, pi/2, 320 );
    
    [~,I] = sort(out2D(:,2),'ascend');
    [~,J1] = sort(out2D(I(1:2),1),'ascend');
    [~,J2] = sort(out2D(I(3:4),1),'descend');
    ID = [I(J1(1)) I(J1(2)) I(J2(1)+2) I(J2(2)+2)];
    rank_xyz = xyz(ID,:);
    rank_rects(vid,:) = reshape(rank_xyz', 1, []);  
end

if length(views)==2
    if any(views==5)
        objects = floorRectangleCombo2Hypothesis( rank_rects, views, room_points );
    else
        objects = rectangleCombo2Hypothesis( rank_rects, views, room_points );
    end
else
    I = find(views~=5);
    objects = rectangleCombo2Hypothesis( rank_rects(I,:), views(I), room_points );
end


end

