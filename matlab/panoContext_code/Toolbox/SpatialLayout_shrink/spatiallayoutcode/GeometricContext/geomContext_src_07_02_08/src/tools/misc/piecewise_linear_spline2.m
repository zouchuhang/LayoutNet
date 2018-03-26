function p = piecewise_linear_spline2(x, y, max_seg)
% performs piecewise linear spline for 2-d points x(1:npts) and
% y(1:ntps) when the number of segments is not known
% p - rows are parameters for each line:
%     y = p(1)*x + p(2) for points with indices p(3) to p(4)     

[x, ind] = sort(x);
y = y(ind);

npts = length(x);

params = get_bottom_hull(x, y);
max_model = min(max_seg, npts-size(params, 1)-1);

mthresh = 0.01;

for m = 1:max_model
    nlines = size(params, 1);    
    gain = zeros(nlines, 1);
    indgain = zeros(nlines, 1);
    splitparams = cell(nlines, 1);
    
    for i = 1:nlines
        % get the gain from splitting a single line
        [gain(i), splitparams{i}, indgain(i)] = get_split_gain(x, y, params(i, :));                        
    end
    [maxgain, maxind] = max(gain);
    %disp(maxgain)
    if indgain(maxind) > mthresh            	
        params = [params(1:maxind-1, :) ; splitparams{maxind} ; params(maxind+1:end, :)];        
    else
        break;
    end
end

p = params;

% if length(x) > 100
% 	figure(1), plot(x, 1-y, '.b');
% 	hold on
% 	for i = 1:size(p, 1)
%         tx(1:2) = x(p(i, 3:4));
%         ty(1:2) = p(i, 1)*tx + p(i, 2);
%         plot(tx, 1-ty, '-g');    
% 	end
% end

% remove edges that are too sharp
p = remove_sharp_edges(p, x, y);
% num_removed = 1;
% while num_removed    
% 	reminds = [];
% 	for i = 1:size(p, 1)
%         if abs(p(i, 1)) > 0.25
%             reminds = [reminds i];
%         end
% 	end
% 	p(reminds, :) = [];
%     num_removed = length(reminds);
% 	for i = 2:size(p, 1)
%         p(i, 3) = p(i-1, 4);
%         i1 = p(i, 3);
%         i2 = p(i, 4);
%         p(i, 1) = (y(i2)-y(i1)) / (x(i2)-x(i1));
%         p(i, 2) = y(i1) - p(i, 1)*x(i1);      
% 	end
% end

p(:, 3) = x(p(:, 3));
p(:, 4) = x(p(:, 4));
    
%disp(p(:, 1)')
if length(x) > 100
	figure(1), hold off, plot(x, 1-y, '.b');
	hold on
	for i = 1:size(p, 1)
        tx(1:2) = p(i, 3:4);
        ty(1:2) = p(i, 1)*tx + p(i, 2);
        plot(tx, 1-ty, '-r');
	end
	hold off
end

%disp(p)

% end piecewise linear regression

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_penalty = get_data_penalty(x, y, params)

data_penalty = 0;
sse = 0;
for i = 1:size(params, 1)
    truey = y(params(i, 3):params(i, 4));
    truex = x(params(i, 3):params(i, 4));
    esty = params(i, 1)*truex + params(i, 2);
    diffy = (truey-esty);
    stdy = std(diffy);
    sse = sse + sum(abs(diffy))/length(y);
end
data_penalty = sse;
%disp(sse)

% end get_data_penalty

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
function [gain, newparams, indgain] = get_split_gain(x, y, params);
% params is params for a single segment
%x = x(params(3):params(4));
%y = y(params(3):params(4));

firstind = params(3);
lastind = params(4);

old_length = length(y);

gain = 0;
indgain = 0;
newparams = params;
x = x(firstind:lastind);
y = y(firstind:lastind);
params(1, 3:4) = [1 length(x)];
currcost = get_data_penalty(x, y, params);
lowcost = currcost;

%step = round(length(x)/100);
step = 1;
for i = 2:step:length(x)-1

    % find best fit and cost from firstind to i and i to lastind 
    tparams(1, 1) = (y(i)-y(1))/(x(i)-x(1));
    tparams(1, 2) = y(1)-x(1)*tparams(1,1);
    tparams(1, 3:4) = [1 i];  
    tparams(2, 1) = (y(length(x))-y(i))/(x(length(x))-x(i));
    tparams(2, 2) = y(i)-x(i)*tparams(2,1);    
    tparams(2, 3:4) = [i length(x)];           
    
    if abs(tparams(1, 1)) < 0.75 & abs(tparams(2, 1)) < 0.75                
        nextcost = get_data_penalty(x, y, tparams);
        %disp([num2str(i) ' ' num2str(nextcost/currcost)])
        
        if nextcost < lowcost
            lowcost = nextcost;  
            newparams = tparams; 
            newparams(:, 3:4) = newparams(:, 3:4) + firstind-1;
        end
    end
end
indgain = currcost - lowcost;
gain = indgain *length(y)/old_length;
% end get_split_gain




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function linp = get_bottom_hull(x, y)
% x is sorted, bottom has higher y
pind = convhull(x, y);
pind = sort(pind(1:end-1));
xhull = x(pind);
yhull = y(pind);

reminds = []; 
for i = 2:length(xhull)-1
    ty = y(1) + (xhull(i)-x(1))*(y(end)-y(1))/(x(end)-x(1));
    if ty > yhull(i)
        reminds = [reminds i];
    end
end
xhull(reminds) = [];
yhull(reminds) = [];
pind(reminds) = [];

linp = zeros(length(xhull)-1, 4);
for i = 1:length(xhull)-1
    linp(i, 1) = (yhull(i+1)-yhull(i)) / (xhull(i+1)-xhull(i));
    if isnan(linp(i, 1))
        disp(['error: ' num2str([xhull(i+1) xhull(i)])])
    end
    linp(i, 2) = yhull(i) - linp(i, 1)*xhull(i); 
    linp(i, 3:4) = [pind(i) pind(i+1)];
end
linp = remove_anomalous_edges(linp, x, y);
linp = remove_sharp_edges(linp, x, y);
% reminds = [];
% for i = 1:size(linp, 1)
%     if abs(linp(i, 1)) > 0.25
%         reminds = [reminds i];
%     end
% end
% linp(reminds, :) = [];
% for i = 2:size(linp, 1)
%     linp(i, 3) = linp(i-1, 4);
%     i1 = linp(i, 3);
%     i2 = linp(i, 4);
%     linp(i, 1) = (y(i2)-y(i1)) / (x(i2)-x(i1));
%     linp(i, 2) = y(i1) - linp(i, 1)*x(i1);     
% end


% if length(x) > 100
% 	figure(1), hold off, plot(x, 1-y, '.b');
% 	hold on
% 	for i = 1:size(linp, 1)
%         tx(1:2) = x(linp(i, 3:4));
%         ty(1:2) = linp(i, 1)*tx + linp(i, 2);
%         plot(tx, 1-ty, '-y');    
% 	end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function linp = remove_anomalous_edges(linp, x, y)

reminds = [];
for i = 2:size(linp,1)-1
    pt = [x(linp(i, 3)) y(linp(i, 3))];
	dists = sqrt((pt(1)-x).^2 + (pt(2)-y).^2);
	if sum(dists<10) <  5
        reminds = [reminds i];
    end
end
linp(reminds, :) = [];
num_removed = length(reminds);
for i = 2:size(linp, 1)
    linp(i, 3) = linp(i-1, 4);
    i1 = linp(i, 3);
    i2 = linp(i, 4);
    linp(i, 1) = (y(i2)-y(i1)) / (x(i2)-x(i1));
    linp(i, 2) = y(i1) - linp(i, 1)*x(i1);      
end            
        
% end remove_anomalous_edges

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = remove_sharp_edges(p, x, y)
% remove edges that are too sharp
num_removed = 1;
while num_removed    
	reminds = [];
	for i = 1:size(p, 1)
        if abs(p(i, 1)) > 0.75
            reminds = [reminds i];
        end
	end
	p(reminds, :) = [];
    num_removed = length(reminds);
	for i = 2:size(p, 1)
        p(i, 3) = p(i-1, 4);
        i1 = p(i, 3);
        i2 = p(i, 4);
        p(i, 1) = (y(i2)-y(i1)) / (x(i2)-x(i1));
        p(i, 2) = y(i1) - p(i, 1)*x(i1);      
	end
end