function [linesout] = LineEndPtsInImg(linesin,w,h)

numlines = size(linesin,1);
linesout = zeros(numlines,4);
x = zeros(numlines,4);
y = zeros(numlines,4);
x(:,1) = 1;
y(:,1) = linesin(:,3) + (x(:,1)-linesin(:,1)).*(linesin(:,4)-linesin(:,3))./(linesin(:,2)-linesin(:,1));

x(:,2) = w;
y(:,2) = linesin(:,3) + (x(:,2)-linesin(:,1)).*(linesin(:,4)-linesin(:,3))./(linesin(:,2)-linesin(:,1));

y(:,3) = 1;
x(:,3) = linesin(:,1) + (y(:,3)-linesin(:,3)).*(linesin(:,2)-linesin(:,1))./(linesin(:,4)-linesin(:,3));

y(:,4) = h;
x(:,4) = linesin(:,1) + (y(:,4)-linesin(:,3)).*(linesin(:,2)-linesin(:,1))./(linesin(:,4)-linesin(:,3));

onbndy = (x>=1 & x<=w & y>=1 & y<=h);
onbndy = onbndy .* repmat([1:4],numlines,1);

[vv,ii] = sort(onbndy,2,'descend');

linesout(:,1) = x(sub2ind(size(x),[1:numlines]',ii(:,1)));
linesout(:,2) = x(sub2ind(size(x),[1:numlines]',ii(:,2)));
linesout(:,3) = y(sub2ind(size(y),[1:numlines]',ii(:,1)));
linesout(:,4) = y(sub2ind(size(y),[1:numlines]',ii(:,2)));


return;
