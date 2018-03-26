function [ obbexp ] = obbExpansion( obb, v )
%OBBEXPANSION Summary of this function goes here
%   Detailed explanation goes here
center = (obb(1:2,:) + obb(5:6,:))/2;
vecx = (obb(3:4,:) + obb(5:6,:))/2 - center;
vecy = (obb(5:6,:) + obb(7:8,:))/2 - center;
disx = sqrt(sum(vecx.^2,1));
disy = sqrt(sum(vecy.^2,1));

unix = vecx./repmat(disx,2,1);
uniy = vecy./repmat(disy,2,1);
for i = 1:size(unix,2)
    if isnan(unix(1,i))
        unix(:,i) = [uniy(2,i);-uniy(1,i)];
    end
    if isnan(uniy(1,i))
        uniy(:,i) = [unix(2,i);-unix(1,i)];
    end
end


disx(disx<0.1) = 0.1;
disx(disx>0.1) = disx(disx>0.1) + v;
disy(disy<0.1) = 0.1;
disy(disy>0.1) = disy(disy>0.1) + v;


% if disx<0.1
%     disx = 0.1;
% else
%     disx = disx + v;
% end
% 
% if disy<0.1
%     disy = 0.1;
% else
%     disy = disy + v;
% end

obbexp = zeros(size(obb));
obbexp(1:2,:) = center - [disx;disx].*unix - [disy;disy].*uniy;
obbexp(3:4,:) = center + [disx;disx].*unix - [disy;disy].*uniy;
obbexp(5:6,:) = center + [disx;disx].*unix + [disy;disy].*uniy;
obbexp(7:8,:) = center - [disx;disx].*unix + [disy;disy].*uniy;

disz = (obb(10,:)-obb(9,:))/2;
centerz = (obb(9,:) + obb(10,:))/2;
disz(disz<0.1) = 0.1;
disz(disz>0.1) = disz(disz>0.1) + v;
obbexp(9,:) = centerz - disz;
obbexp(10,:) = centerz + disz;

end

