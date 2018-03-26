function rinds = APPspInds2RegionInds(map, sinds)
% rinds = APPspInds2RegionInds(map, sinds)
% Gets the sp in each region
%
% Copyright(C) Derek Hoiem, Carnegie Mellon University, 2005
% Permission granted to non-commercial enterprises for
% modification/redistribution under GNU GPL.  
% Current Version: 1.0  09/30/2005

nr = max(map);
rinds = cell(nr, 1);
for r = 1:nr
    count = 0;
    rs = find(map==r);
    for k = 1:length(rs)
        count = count + length(sinds{rs(k)});
    end         
    rinds{r} = zeros(count, 1);
    
    count = 0;
    for k = 1:length(rs) 
        rinds{r}(count+1:count+length(sinds{rs(k)})) = sinds{rs(k)};
        count = count + length(sinds{rs(k)});
    end
end