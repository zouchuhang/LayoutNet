function [ indicator ] = insideIndicator( hyps, coor )
%INSIDEINDICATOR Summary of this function goes here
%   Detailed explanation goes here
num_coor = size(coor,1);
num_object = size(hyps,1);
indicator = false(num_coor, num_object);

for i = 1:num_object
    if hyps(i).type==1
        p = hyps(i).out_points_w([1 2 3 4],:);
    else
        p = hyps(i).align;
    end
    p = p./repmat(sqrt(sum(p.^2,2)),1,3);
    vp = sum(p,1); vp = vp./norm(vp); %vp(3) = 0; 
    [out2D, valid, ~] = projectPoint2SeparateView( p, vp, pi/3, 100 );
    if any(~valid)
        vp = sum(p,1); vp(1:2) = 0; vp = vp./norm(vp); 
        [out2D, valid, ~] = projectPoint2SeparateView( p, vp, pi/3, 100 );
        if any(~valid)
            vp = sum(p,1); vp(3) = 0; vp = vp./norm(vp); 
            [out2D, valid, ~] = projectPoint2SeparateView( p, vp, pi/3, 100 );
        end
    end
    K = convhull( out2D(:,1), out2D(:,2));
    [ inside, ~, ~ ] = insideCone( p(K(end:-1:2),:), coor, 0 );
    indicator(:,i) = inside;
end

end

