function [colinear]=chckcollinearity(points,thres)

if ~exist('thres','var')
   thres =10^(-5);
end

x=points(:,1);
y=points(:,2);


if (numel(x)*sum(x.*x) == (sum(x))^2)
    colinear=1;       %vertical line
    return;
end
% fit line by lls
Mat1=inv([numel(x) sum(x); sum(x) sum(x.*x)]);
vec1=[sum(y); sum(x.*y)];

linepar=Mat1*vec1;


fiterr=sum((y-(linepar(1)+linepar(2)*x)).^2);

if fiterr < thres
    colinear=1;
else
    colinear=0;
end
return;