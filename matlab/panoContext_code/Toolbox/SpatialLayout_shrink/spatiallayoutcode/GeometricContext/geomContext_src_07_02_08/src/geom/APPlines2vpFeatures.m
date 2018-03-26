function data = APPlines2vpFeatures(lines, max_angle)
% APPlines2vpFeatures(lines, max_angle)
% Computes vanishing point related features from straight lines
% data is currently 19 elements, line histograms not used
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005

nbins=0;
data = zeros(nbins+19, 1);

nLines = size(lines, 1);

if nLines <= 1
    data(nbins+(1:2)) = 0.5;
    data(nbins+3) = 0;
    data(nbins +(4:10)) = 0;
    return;
end

a = lines(:, 5);
x1 = lines(:, 1);
y1 = lines(:, 3);
x2 = lines(:, 2);
y2 = lines(:, 4);
r = lines(:, 6);

slength = sqrt((x2-x1).^2 + (y2-y1).^2);

xc = (x1+x2)/2;
yc = (y1+y2)/2;

m = (y2-y1)./(x2-x1+1E-10);

a1 = atan(tan(a));


adist = pdist(a, 'cityblock');
%adist = mod(adist, pi/2) < max_angle;
adist = adist < max_angle;
adistM = squareform(adist);


% get |m1-m2|/|x1-x2|, where the direction of x is perpendicular to m 
parabins = zeros(8, 1);
paracount = zeros(8, 1);
parabins2 = zeros(8, 1);
leftcount = 0.05;
rightcount = 0.05;
upcount = 0.05;
downcount = 0.05;
histx = [(-pi/2):(pi/8):(pi/2)]-pi/16;
for i = 1:nLines    
    bn = sum(a1(i) > histx);
    if bn==9
        bn = 1;
    end

    m1 = m(i);  x1 = xc(i);  y1 = yc(i); 
    %disp(num2str([x1 y1]))
    isvert = abs(sin(a(i))) > abs(cos(a(i)));
    ishorz = abs(cos(a(i))) > abs(sin(a(i)));
    for j = i+1:nLines
        if adistM(i, j) 
            sl = min(slength([i j]));
            
            m2 = m(j); x2 = xc(j);  y2 = yc(j);
           
            if ~(m1==m2 & (y2-y1)==m1*(x2-x1)) % check for colinearity
                
                xi = (y2-y1+m1*x1-m2*x2)/(m1-m2+1E-10);
                yi = (-y1*m2+m1*y2-m1*m2*x2+m1*x1*m2)/(m1-m2+1E-10);                              
                
                if sqrt(xi.^2+yi.^2) > 1 | m1==m2
                    parabins(bn) = parabins(bn) + sl;
                end
                paracount(bn) = paracount(bn) + sl;    
                
                if sqrt(xi.^2+yi.^2) > 3.5 | m1==m2
                    parabins2(bn) = parabins2(bn) + sl;
                end
                
                if isvert
                    if yi > (y1+y2)/2
                        upcount = upcount + sl;
                    else
                        downcount = downcount + sl;
                    end
                else
                    if xi > (x1+x2)/2
                        rightcount = rightcount + sl;
                    else
                        leftcount = leftcount + sl;
                    end
                end
            end
        end
    end
end

ind = find(paracount>0);
parabins(ind) = parabins(ind) ./ paracount(ind);
parabins2(ind) = parabins2(ind) ./ paracount(ind);
%data(1:nbins) = angle_hist;  % extent of different vanishing point directions
%data(nbins+1) = angle_entropy; % measure of how many different vanishing points there may be
data(1) = rightcount ./ (leftcount+rightcount); % differentiate between face right and face left
data(2) = upcount ./ (upcount+downcount); % differentiate between face up and face down
%data(3) = sum(slength); % total number of line pixels in region
data(3) = mean(adist);  % percent of pairs of angles that are roughly parallel
data(4:11) = parabins; % measure of "parallelness" at different angles
data(12:19) = parabins2;

%disp(['bins: ' num2str(data(1:nbins)')]);               
%disp(['rand: ' num2str(data(nbins+1:nbins+5)')])
%disp(['para: ' num2str(data(nbins+6:end)')])
            
            
            


