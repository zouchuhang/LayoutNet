function p = piecewise_linear_spline(x, y, max_seg)
% performs piecewise linear spline for 2-d points x(1:npts) and
% y(1:ntps) when the number of segments is not known
% p - rows are parameters for each line:
%     y = p(1)*x + p(2) for points with indices p(3) to p(4)     

[x, ind] = sort(x);
y = y(ind);

npts = length(x);
max_model = min(max_seg, npts-1);
mdl_cost = zeros(max_model, 1);

if npts==1
    p = [0 y x x];
end

% y = p(1)*x + p(2)
linp(1, 1) = (y(npts) - y(1)) / (x(npts) - x(1)); 
linp(1, 2) = y(1) - x(1)*linp(1,1); 
linp(1, 3:4) = [1 npts]; % endpoint indices
params = linp;

mthresh = 0.005*max(y);

m = 1;

while 1
    
    m = m + 1;
    
    gain = zeros(m-1, 1);
    indgain = zeros(m-1, 1);
    splitparams = cell(m-1, 1);
    
    for i = 1:size(params, 1)
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
    if m == max_model
        break;
    end
end

p = params;
p(:, 3) = x(p(:, 3));
p(:, 4) = x(p(:, 4));

% if length(x) > 100
% figure(1), hold off, plot(x, 1-y, '.b');
% hold on
% for i = 1:size(p, 1)
%     tx(1:2) = p(i, 3:4);
%     ty(1:2) = p(i, 1)*tx + p(i, 2);
%     plot(tx, 1-ty, '-r');
% end
% hold off
% end



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








